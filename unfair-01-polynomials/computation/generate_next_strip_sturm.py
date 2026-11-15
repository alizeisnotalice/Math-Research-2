#!/usr/bin/env python3
"""Exact finite-r ledger for m=2r+2 endpoint prefixes.

For r>=9 feasibility forces 13/20<s<17/20 (proved separately).  After
b_(r+2)=0, t=U_(r+2)/D.  The only c-prefix branch is u=c_(2r).
u=1 requires F_r=b_(2r+1)=0; u=0 requires F_r+s=0.

This script uses SymPy's exact QQ Sturm root counter only.  It writes an
atomic JSON checkpoint after every r.  A zero count is a rigorous exclusion
of that branch conditional on the proved interval lemma; a nonzero count is
only a survivor to be tested with the remaining inequalities/later prefix.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
import sympy as S

s = S.symbols("s")
LO, HI = S.Rational(13, 20), S.Rational(17, 20)
OUT = Path("/Volumes/data/CodexTemp/T2/unfair-01-polynomials/round011/generated-endpoint-ledger.json")


def atomic_write(rows):
    assert Path("/Volumes/data").is_mount()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    target = OUT.with_suffix(".json.tmp")
    target.write_text(json.dumps(rows, ensure_ascii=False, indent=2) + "\n")
    os.replace(target, OUT)


def poly_num(expr):
    return S.Poly(S.together(expr).as_numer_denom()[0], s, domain=S.QQ)


def row_for(r: int):
    D = 1 - 2*s + 3*s*s
    U = lambda n: sum((-s)**j for j in range(n+1))
    W = lambda n: sum((j+1)*(-s)**j for j in range(n+1))
    t = S.cancel(U(r+2)/D)
    # all-one prefix through 2r+1; it is the u=1 branch.
    F = S.cancel(U(2*r+1)-t*W(r+1)+t*t*(1-3*s))
    ans = {"r": r, "interval": [str(LO), str(HI)], "branches": {}}
    for branch, expr in (("u1", F), ("u0", S.cancel(F+s))):
        p = poly_num(expr)
        # SymPy's endpoint convention is immaterial only after this exact check.
        assert p.eval(LO) != 0 and p.eval(HI) != 0
        count = int(p.count_roots(LO, HI))
        entry = {
            "degree": p.degree(),
            "sturm_roots_open_interval": count,
            "numerator": str(p.as_expr()),
        }
        if count:
            # disjoint rational isolators make the outstanding roots explicit.
            intervals = []
            for (a, b), mult in p.intervals(eps=S.Rational(1, 10)**18):
                if LO < a and b < HI:
                    intervals.append({"lo": str(a), "hi": str(b), "mult": int(mult)})
            entry["isolators"] = intervals
        ans["branches"][branch] = entry
    return ans


def main():
    rows = []
    if OUT.exists():
        rows = json.loads(OUT.read_text())
    done = {row["r"] for row in rows}
    for r in range(9, 72):
        if r in done:
            continue
        record = row_for(r)
        rows.append(record)
        rows.sort(key=lambda x: x["r"])
        atomic_write(rows)
        print(json.dumps(record, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
