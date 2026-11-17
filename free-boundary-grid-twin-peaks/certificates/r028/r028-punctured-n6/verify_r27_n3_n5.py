#!/usr/bin/env python3
"""Cross-check the n=6 C++ kernel against R027's Python exact tables."""
import json
import subprocess
import sys
from pathlib import Path


def main(binary: Path, r027: Path) -> None:
    for n in (3, 4, 5):
        old = json.loads((r027 / f"punctured-table-n{n}.json").read_text())
        new = json.loads(subprocess.check_output([str(binary), str(n)], text=True))
        assert int(new["A2"]) == old["A2_total"]
        assert int(new["B"]) == old["B_total"]
        assert [int(x["A2_v"]) for x in new["per_vertex"]] == [x["A2_v"] for x in old["per_vertex"]]
        assert [int(x["B_v"]) for x in new["per_vertex"]] == [x["B_v"] for x in old["per_vertex"]]
        print(f"n={n}: A2, B, and every vertex matched R027")


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))
