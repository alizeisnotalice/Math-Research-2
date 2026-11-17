#!/usr/bin/env python3
"""Finite exact checks of the proposed quadratic tail envelope, n=3..5."""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path


def from_full(path: Path):
    x = json.loads(path.read_text())
    h = [[int(a) for a in q["histogram_L"]] for q in x["holes"]]
    B = sum(int(q["B_v_from_D"]) for q in x["holes"])
    T = [sum(sum(a[t - 1:]) for a in h) for t in range(1, len(h[0]) + 1)]
    return B, T


def from_n5(path: Path):
    x = json.loads(path.read_text())
    return int(x["B_total"]), list(map(int, x["tail_T_ge_t"]))


def measure(n: int, B: int, tails: list[int]):
    vals = [(Fraction(t * (t + 1) * q, B), t) for t, q in enumerate(tails, 1) if q]
    q, t = max(vals)
    return {"n": n, "B": B, "max_t_times_t_plus_1_T_over_B": str(q),
            "max_decimal": float(q), "attained_t": t,
            "quadratic_envelope_constant_8_holds": q <= 8}


def main(root):
    rows = [measure(3, *from_full(root / "tail-full-n3.json")),
            measure(4, *from_full(root / "tail-full-n4.json")),
            measure(5, *from_n5(root / "n5-tail-summary.json"))]
    assert all(r["quadratic_envelope_constant_8_holds"] for r in rows)
    out = {"status": "passed", "scope": "Exact finite n=3..5 checks only; not a general bound.",
           "candidate": "T_t <= 8 B/[t(t+1)] for t>=2", "rows": rows}
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    import argparse
    p=argparse.ArgumentParser(); p.add_argument("--r29",type=Path,required=True)
    main(p.parse_args().r29)
