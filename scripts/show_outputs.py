"""Print each code cell's console output, or its error, for a notebook."""
import json, sys

def text(chunk):
    return "".join(chunk) if isinstance(chunk, list) else chunk

for path in sys.argv[1:]:
    print(f"===== {path}")
    nb = json.load(open(path))
    for i, cell in enumerate(nb["cells"]):
        if cell["cell_type"] != "code":
            continue
        errors = [o for o in cell.get("outputs", []) if o.get("output_type") == "error"]
        if errors:
            e = errors[0]
            print(f"[{i}] ERROR {e['ename']}: {text(e['evalue'])[:200]}")
            continue
        out = "".join(
            text(o.get("text", ""))
            for o in cell.get("outputs", [])
            if o.get("output_type") == "stream"
        ).strip()
        if out:
            print(f"[{i}] " + out.replace("\n", "\n    "))
