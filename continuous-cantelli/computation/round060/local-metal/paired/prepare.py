"""Take the stored representative indices from both Metal scans and the old audit."""
import json
from pathlib import Path
from sys import argv


def read(path):
    return json.loads(path.read_text())


def main():
    old = Path(argv[1])
    new = Path(argv[2])
    chosen = {}

    def append(origin, row):
        ident = int(row["index"])
        selected = chosen.setdefault(ident, {"index": ident, "point": row["point"], "origin": origin})
        if selected["origin"] != origin:
            selected["origin"] += "," + origin

    for label, scan in (("128 new",new),("96 old",old)):
        for name, row in read(scan / "quartic-main.json")["summaries"].items():
            append(f"{label} {name}", row["best"])
    audit = read(old / "audit/audit-results.json")
    append("previously unstable",audit["previously_unstable"])
    for name, modes in audit["strata"].items():
        for mode, cases in modes.items():
            for row in cases:
                append(f"old sparse {name} {mode}",row)
    for row in read(old / "rare/rare-results.json")["samples"]:
        append("old rare flip",row)
    points = sorted(chosen.values(), key=lambda x: x["index"])
    Path(argv[3]).write_text(json.dumps(points,indent=2,sort_keys=True)+"\n")


if __name__ == "__main__":
    main()
