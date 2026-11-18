"""Replay the exact convolution verifier on all archived Metal results."""

import argparse
from fractions import Fraction
import json
from pathlib import Path
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--scratch", required=True, type=Path)
    args = parser.parse_args()
    verifier = Path(__file__).resolve().parent.parent / "063-gap7-numeric" / "verify_metal_samples.py"
    args.scratch.mkdir(parents=True, exist_ok=True)
    samples = []
    for source in sorted(args.results.glob("m*n*.json")):
        data = json.loads(source.read_bytes())
        assert data["status"] == "finished_numeric_search"
        m, n = data["m"], data["n"]
        assert source.name == f"m{m}n{n}.json"
        assert data["trajectories"] == 32768 and data["completed"] == data["steps"] == 32000
        out = args.scratch / source.name
        subprocess.run(
            [sys.executable, "-B", str(verifier), str(source), "--out", str(out)],
            check=True, stdout=subprocess.DEVNULL,
        )
        report = json.loads(out.read_bytes())
        assert len(report["residuals"]) == 16
        assert all(Fraction(value) > 0 for value in report["residuals"])
        samples.append(dict(m=m, n=n, trajectories=data["trajectories"], steps=data["steps"], seed=data["seed"], **report))
    assert len(samples) == 40 and len({row["seed"] for row in samples}) == 40
    report = dict(status="verified_rounded_samples", pairs=len(samples), samples=sum(row["samples"] for row in samples), rows=samples)
    with args.out.open("x") as output:
        json.dump(report, output, indent=2)
        output.write("\n")
    print(report["pairs"], report["samples"])


if __name__ == "__main__":
    main()
