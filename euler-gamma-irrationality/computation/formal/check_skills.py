"""Cross-check the user's independently installed gamma skill entry points."""

import argparse
import json
from fractions import Fraction
from pathlib import Path
from subprocess import check_output
import sys


def run(root, route):
    name, script, options = {
        "apt": (
            "math-gamma-multiple-laguerre-approximants",
            "comparison_pairs.py",
            ["aptekarev"],
        ),
        "mixed": ("math-gamma-multiple-laguerre-approximants", "laguerre_pairs.py", []),
        "riv": (
            "math-gamma-multiple-laguerre-approximants",
            "comparison_pairs.py",
            ["rivoal"],
        ),
        "cmf": ("math-gamma-meijer-cmf", "cmf_pairs.py", []),
    }[route]
    rows = json.loads(Path(__file__).with_name(route + ".json").read_text())
    if route == "cmf":
        rows = rows[0]["pairs"]
    other = json.loads(
        check_output(
            [
                sys.executable,
                "-B",
                str(root / name / "scripts" / script),
                *options,
                *[str(r["n"]) for r in rows],
            ]
        )
    )
    assert len(other) == len(rows)
    for a, b in zip(rows, other):
        assert a["n"] == b["n"]
        if route in ("mixed", "cmf"):
            assert (a["P"], a["Q"]) == (b["P"], b["Q"])
        else:
            assert Fraction(a["P"], a["Q"]) == Fraction(b["p"]) / Fraction(b["q"])
            if route == "apt":
                assert (a["P"], a["Q"]) == (Fraction(b["p"]), Fraction(b["q"]))
    print(route, len(rows), "OK")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("route", nargs="+", choices=("apt", "mixed", "riv", "cmf"))
    ap.add_argument(
        "--skills", type=Path, default=Path("/Volumes/data/CodexData/skills")
    )
    args = ap.parse_args()
    for route in args.route:
        run(args.skills, route)
