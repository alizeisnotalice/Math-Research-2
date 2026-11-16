#!/usr/bin/env python3
"""Strict Decimal audit for r021 higher inverse-fibre SMC CSVs.

The reported delta quantities use independent replicate residuals. They are
diagnostics for this finite randomized computation, not confidence guarantees.
"""
import argparse
import csv
import hashlib
import json
from decimal import Decimal as D, localcontext
from pathlib import Path


def summarize(path):
    meta = json.loads(Path(str(path) + ".metadata.json").read_text())
    with path.open() as f:
        rows = list(csv.DictReader(f))
    n, P, R, L, T = (int(meta[k]) for k in ("n", "P", "R", "L", "threads"))
    assert n >= 2 and P > 0 and R > 1 and L >= 0 and 1 <= T <= 12
    assert len(rows) == R and {int(r["rep"]) for r in rows} == set(range(R))
    tau, seedbase = n * n + 1, int(meta["seed_base"])
    with localcontext() as ctx:
        ctx.prec = 90
        def close(x, y):
            assert abs(x-y) <= D("5e-14") * max(abs(x), abs(y)), (x, y)
        for r in rows:
            rep = int(r["rep"])
            assert (int(r["n"]), int(r["tau"]), int(r["L"])) == (n, tau, L)
            assert int(r["seed"]) == (seedbase + tau*1000003 + rep*1009 + L*1000000007) % (1 << 64)
            assert 0 <= int(r["accept"]) <= int(r["proposals"])
            assert int(r["proposals"]) == P * L * (n*n-1)
            assert int(r["dead"]) >= 0
            z = D(r["Z"]); assert z.is_finite() and z >= 0
            for d in range(2, 6):
                m = D(r[f"Zmean_m{d}"])
                assert m.is_finite() and m >= 0
                if z: close(m/z, D(r[f"weighted_ratio_m{d}"]))
                else: assert m == 0
        zs = [D(r["Z"]) for r in rows]
        total = sum(zs); assert total > 0
        mean = total / R
        result = {
            "path": path.name,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "metadata": meta,
            "replicates": R,
            "mean_Z": str(mean),
            "mean_Z_delta_SE": str((sum((z-mean)**2 for z in zs) / D(R*(R-1))).sqrt()),
            "replicate_weight_ESS": str(total**2 / sum(z*z for z in zs)),
            "largest_replicate_weight_share": str(max(zs) / total),
            "zero_normalizer_replicates": sum(z == 0 for z in zs),
            "dead_particle_rank_events": sum(int(r["dead"]) for r in rows),
            "minimum_particle_parent_ESS": min(float(r["min_parent_ess"]) for r in rows),
            "proposals": sum(int(r["proposals"]) for r in rows),
            "accepted": sum(int(r["accept"]) for r in rows),
            "moments": {},
            "scope": "Finite SMC diagnostic only; no certified confidence interval, mixing bound, or asymptotic conclusion.",
        }
        for d in range(2, 6):
            ms = [D(r[f"Zmean_m{d}"]) for r in rows]
            ratio = sum(ms) / total
            residual = [m - ratio*z for m, z in zip(ms, zs)]
            se = (sum(x*x for x in residual) / D(R*(R-1))).sqrt() / mean
            scale = D(n*n) ** (d-1)
            result["moments"][str(d)] = {
                "ratio_of_means": str(ratio), "delta_SE": str(se),
                "scaled_by_N_to_d_minus_1": str(ratio / scale),
                "scaled_delta_SE": str(se / scale),
            }
        return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("csv", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = summarize(args.csv)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
