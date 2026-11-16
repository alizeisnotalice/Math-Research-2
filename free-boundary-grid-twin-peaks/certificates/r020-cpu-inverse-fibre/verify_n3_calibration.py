#!/usr/bin/env python3
"""Strict fibre-table crosscheck plus deterministic n=3 SMC calibration gate."""
import csv
from fractions import Fraction
import json
from pathlib import Path


BASE = Path(__file__).parent
A1, A2D2, A2D3 = 29568, 63552, 64512
TARGET2, TARGET3 = Fraction(A2D2, A1), Fraction(A2D3, A1)


def csv_rows(name):
    with (BASE / name).open() as f:
        return list(csv.DictReader(f))


def main():
    # The two all-word tables are independently produced and compared exactly.
    assert (BASE / "n3-python-fibres.csv").read_bytes() == (BASE / "n3-cpp-fibres.csv").read_bytes()
    table = csv_rows("n3-cpp-fibres.csv")
    assert len(table) == A1
    assert sum(int(r["m2"]) for r in table) == A2D2
    assert sum(int(r["m3"]) for r in table) == A2D3
    runs = []
    # These finite-MC gates are intentionally loose diagnostics, not confidence
    # statements: the exact all-word table above is the actual correctness test.
    for L in (0, 4, 16):
        rows = csv_rows(f"n3-p17-r512-l{L}.csv")
        assert len(rows) == 512
        assert all(int(r["tau"]) == 10 and int(r["dead"]) == 0 for r in rows)
        z = sum(float(r["Z"]) for r in rows)
        zm2 = sum(float(r["Zmean_m2"]) for r in rows)
        zm3 = sum(float(r["Zmean_m3"]) for r in rows)
        mean_z, ratio2, ratio3 = z / len(rows), zm2 / z, zm3 / z
        errors = (abs(mean_z/A1 - 1), abs(ratio2/float(TARGET2) - 1),
                  abs(ratio3/float(TARGET3) - 1))
        assert max(errors) < .10, (L, errors)
        runs.append({"L": L, "replicates": len(rows), "mean_Z": mean_z,
                     "weighted_ratio_m2": ratio2, "weighted_ratio_m3": ratio3,
                     "relative_errors_vs_exact": list(errors)})
    out = {"status": "passed", "exact_table": {"A1": A1, "sum_m2": A2D2,
           "sum_m3": A2D3, "target_m2": str(TARGET2), "target_m3": str(TARGET3)},
           "smc_gate": {"replicates": 512, "relative_tolerance": .10,
                        "scope": "deterministic finite-MC diagnostic only"}, "runs": runs}
    (BASE / "n3-calibration-verified.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
