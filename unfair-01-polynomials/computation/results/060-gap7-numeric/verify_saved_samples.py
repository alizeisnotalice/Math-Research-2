"""Recompute residuals for all retained rounded GPU samples."""

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path


def scaled(coefficients):
    ratios = [value.as_integer_ratio() for value in coefficients]
    denominator = max(divisor for _, divisor in ratios)
    integers = [value * (denominator // divisor) for value, divisor in ratios]
    assert all(denominator % divisor == 0 for _, divisor in ratios)
    assert integers[0] == integers[-1] == denominator
    assert all(0 <= value <= denominator for value in integers)
    return integers, denominator


def residual(case, m, n):
    a, da = scaled(case["A"])
    b, db = scaled(case["B"])
    assert (len(a), len(b)) == (m + 1, n + 1)
    denominator = da * db
    coefficients = [0] * (m + n + 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            coefficients[i + j] += ai * bj
    nearest = [int(2 * coefficient >= denominator) for coefficient in coefficients]
    assert nearest[:9] == [1] * 7 + [0, 1]
    error = Fraction(max(abs(x - y * denominator) for x, y in zip(coefficients, nearest)), denominator)
    assert error > 0
    assert abs(float(error) - case["residual_linf"]) <= 1e-12
    return error


def audit(directory):
    source = (directory / "result.json").read_bytes()
    result = json.loads(source)
    inputs = json.loads((directory / "input.json").read_bytes())["parameters"]
    release = json.loads((directory / "release.json").read_bytes())
    assert result["status"] == release["terminal_status"] == "finished_numeric_search"
    assert result["steps"] == release["step"] == inputs["steps"]
    assert release["cleanup_called"]
    windows = [window["cases"] for window in result["completed_windows"]] + [result["current"]]
    assert len(windows) == inputs["steps"] // inputs["restart"]
    assert all(len(window) == 16 for window in windows)
    errors = [residual(case, inputs["m"], inputs["n"]) for window in windows for case in window]
    assert all(point["device"] == "mps:0" for point in result["history"])
    return {"path": directory.name, "sha256": hashlib.sha256(source).hexdigest(),
            "retained_samples": len(errors), "exact_minimum_recorded_residual": str(min(errors)),
            "minimum_recorded_residual": float(min(errors)),
            "seconds": result["history"][-1]["seconds"], "cleanup_pid": release["pid"]}


def main():
    assert __debug__
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("directories", nargs="+", type=Path)
    args = parser.parse_args()
    report = {"status": "verified_rounded_samples", "cases": [audit(path) for path in args.directories]}
    with args.out.open("x") as output:
        json.dump(report, output, indent=2)
        output.write("\n")
    print(json.dumps(report))


if __name__ == "__main__":
    main()
