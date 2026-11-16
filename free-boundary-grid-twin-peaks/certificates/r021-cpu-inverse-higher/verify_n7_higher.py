#!/usr/bin/env python3
"""Audit one r021 n=7 SMC CSV against the independently exact r018 histogram.

This checks file completeness and deterministic metadata, then forms the
ratio-of-means estimator sum_r Z_r m_{d,r}/sum_r Z_r.  It deliberately does
not turn agreement of this finite Monte Carlo run into an asymptotic claim.
"""
import csv
import json
import math
import sys
from decimal import Decimal, getcontext
from pathlib import Path

getcontext().prec = 80
HERE = Path(__file__).resolve().parent
CSV = HERE / "n7-p4096-r32-l16.csv"
META = Path(str(CSV) + ".metadata.json")
EXACT = HERE / "reference-r018-n7-distance-histogram.json"
OUT = HERE / "n7-higher-exact-crosscheck.json"

def fail(msg):
    raise SystemExit("FAILED: " + msg)

def main():
    meta = json.loads(META.read_text())
    exact = json.loads(EXACT.read_text())
    with CSV.open(newline="") as f:
        rows = list(csv.DictReader(f))
    required = {"n", "P", "R", "L", "threads", "seed_base", "tau"}
    if required - set(meta):
        fail("metadata keys missing")
    n, P, R, L, T, tau = (int(meta[k]) for k in ("n", "P", "R", "L", "threads", "tau"))
    base = int(meta["seed_base"])
    if (n, P, R, L, T, tau) != (7, 4096, 32, 16, 12, 50):
        fail("unexpected run parameters")
    if len(rows) != R:
        fail(f"need {R} rows, found {len(rows)}")
    got_reps = {int(r["rep"]) for r in rows}
    if got_reps != set(range(R)):
        fail("replicate ids are incomplete or duplicated")
    expected_props = P * L * (n * n - 1)
    for r in rows:
        rep = int(r["rep"])
        expected_seed = (base + tau * 1000003 + rep * 1009 + L * 1000000007) % (1 << 64)
        if int(r["seed"]) != expected_seed:
            fail(f"bad seed at rep {rep}")
        if int(r["dead"]) != 0 or int(r["proposals"]) != expected_props:
            fail(f"bad dead/proposals record at rep {rep}")
        z = float(r["Z"])
        if not (math.isfinite(z) and z > 0):
            fail(f"nonpositive/nonfinite Z at rep {rep}")
        for d in range(2, 6):
            q = float(r[f"weighted_ratio_m{d}"])
            direct = float(r[f"Zmean_m{d}"]) / z
            if not math.isclose(q, direct, rel_tol=3e-13, abs_tol=0.0):
                fail(f"stored ratio mismatch for d={d}, rep={rep}")
    den = sum(Decimal(r["Z"]) for r in rows)
    result = {
        "scope": "finite n=7 SMC calibration only; no mixing or asymptotic conclusion",
        "csv": CSV.name, "csv_rows": len(rows), "metadata": meta,
        "exact_histogram_source": EXACT.name,
        "checks": {"complete_replicates": True, "seed_formula": True,
                   "dead_zero": True, "proposals_per_rep": expected_props,
                   "stored_ratio_identity": True},
        "ratio_of_means": {},
    }
    for d in range(2, 6):
        estimate = sum(Decimal(r[f"Zmean_m{d}"]) for r in rows) / den
        target = Decimal(exact["distance_counts"][str(d)]) / Decimal(exact["A1"])
        rel = estimate / target - 1
        result["ratio_of_means"][str(d)] = {
            "estimate": str(estimate), "exact_A2_distance_over_A1": str(target),
            "relative_error": str(rel),
        }
    OUT.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))

if __name__ == "__main__":
    main()
