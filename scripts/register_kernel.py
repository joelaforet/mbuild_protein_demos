"""Register a Jupyter kernel for this project's pixi environment.

Run it through pixi so it sees the environment it should register::

    pixi run kernel              # the default environment
    pixi run -e cuda13 kernel    # a GPU environment

It writes one kernelspec into the current user's Jupyter data directory
(``~/.local/share/jupyter`` on Linux, ``~/Library/Jupyter`` on macOS), so
nothing inside the project changes. The kernel's command is ``pixi run
... python -m ipykernel_launcher``: starting through ``pixi run``
activates the environment and returns, so the kernel has the
environment's PATH and ``packmol`` and the OpenMM plugins are found.

The kernelspec also names the environment's interpreter in its metadata.
The VS Code Jupyter extension uses that to find nglview's widget
JavaScript in the environment; without it nglview reports "No version of
module nglview-js-widgets is registered". Naming the interpreter makes VS
Code ask the Python extension to activate the environment, and that step
hangs on ``pixi shell`` unless ``scripts/pixi-vscode`` is in place (see
``.vscode/settings.json`` and the README's VS Code section).

``--remove`` deletes the kernelspec.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
from pathlib import Path

from jupyter_core.paths import jupyter_data_dir

PROJECT = Path(__file__).resolve().parent.parent
MANIFEST = PROJECT / "pixi.toml"
BASE_NAME = "mbuild-protein-demos"


def kernel_identity() -> tuple[str, str]:
    """Return (kernelspec name, display name) for the active pixi env."""
    env = os.environ.get("PIXI_ENVIRONMENT_NAME", "default")
    if env == "default":
        return BASE_NAME, "Python (mbuild_protein_demos)"
    return f"{BASE_NAME}-{env}", f"Python (mbuild_protein_demos {env})"


def find_pixi() -> str:
    exe = os.environ.get("PIXI_EXE") or shutil.which("pixi")
    if not exe:
        sys.exit("pixi not found. Run this as `pixi run kernel`.")
    return exe


def write_kernelspec(name: str, display: str, env: str) -> Path:
    kernel_dir = Path(jupyter_data_dir()) / "kernels" / name
    kernel_dir.mkdir(parents=True, exist_ok=True)
    argv = [
        find_pixi(), "run", "--quiet",
        "--manifest-path", str(MANIFEST),
        "--environment", env,
        "python", "-Xfrozen_modules=off", "-m", "ipykernel_launcher",
        "-f", "{connection_file}",
    ]
    spec = {
        "argv": argv,
        "display_name": display,
        "language": "python",
        "interrupt_mode": "signal",
        "metadata": {
            "debugger": True,
            # Lets the VS Code Jupyter extension locate the environment's
            # widget JavaScript (nglview) via the interpreter's prefix.
            "interpreter": {"path": sys.executable},
        },
    }
    (kernel_dir / "kernel.json").write_text(json.dumps(spec, indent=1) + "\n")
    return kernel_dir


def remove(name: str) -> None:
    kernel_dir = Path(jupyter_data_dir()) / "kernels" / name
    if kernel_dir.exists():
        shutil.rmtree(kernel_dir)
        print(f"removed {kernel_dir}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--remove", action="store_true",
                        help="unregister the kernel")
    args = parser.parse_args()

    if not MANIFEST.exists():
        sys.exit(f"no pixi.toml at {MANIFEST}")
    env = os.environ.get("PIXI_ENVIRONMENT_NAME", "default")
    name, display = kernel_identity()

    if args.remove:
        remove(name)
        return

    kernel_dir = write_kernelspec(name, display, env)
    print(f"kernel '{display}' registered at {kernel_dir}")
    print(
        "\nIn VS Code: open this folder as the workspace, reload the window "
        "once (Developer: Reload Window),\nthen in a notebook choose\n"
        f"  Select Kernel -> Jupyter Kernel... -> {display}"
    )


if __name__ == "__main__":
    main()
