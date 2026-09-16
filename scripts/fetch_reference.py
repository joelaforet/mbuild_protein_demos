"""Download the workshop's prepared 7KI0 structure.

The file is 21 MB and is not committed. It is only needed to regenerate
``semaglutide_reference.pdb`` and ``semaglutide_apo.pdb``, which are.
The notebooks never read it.
"""

import urllib.request
from pathlib import Path

URL = (
    "https://raw.githubusercontent.com/openforcefield/2026-virtual-workshops"
    "/main/ptm/7KI0_prepared.pdb"
)
TARGET = Path(__file__).resolve().parent.parent / "assets_cache" / "7KI0_prepared.pdb"


def fetch():
    if TARGET.exists():
        print(f"{TARGET.name} is already here")
        return
    TARGET.parent.mkdir(parents=True, exist_ok=True)
    print(f"downloading {URL}")
    urllib.request.urlretrieve(URL, TARGET)
    print(f"wrote {TARGET} ({TARGET.stat().st_size // 1_000_000} MB)")


if __name__ == "__main__":
    fetch()
