#!/Users/zzh/.venvs/main/bin/python
"""Exact conditional ledger for the degree-six linked {2,3} core.

This script does not prove the root-envelope bounds used below.  It merely
enumerates exactly the finite parameter box conditional on those bounds, and
does not enumerate the 645 possible (d,e) pairs for each surviving row.
"""
from __future__ import annotations

import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from math import isqrt
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent


def nonsquare(n: int) -> bool:
    return isqrt(n) ** 2 != n


def strict_two_embedding_bound(a0: int, a1: int, delta: int, bound: int) -> bool:
    """Return max(|a0+a1 sqrt(delta)|,|a0-a1 sqrt(delta)|)<bound exactly."""
    rem = bound - abs(a0)
    return rem > 0 and a1 * a1 * delta < rem * rem


def main() -> None:
    stages: Counter[str] = Counter()
    rows: list[dict[str, int]] = []

    # Conditional source bounds supplied by the parent task:
    # b != 0, |b|<=7; |T-b^2|+sqrt(delta)<51/4;
    # 8|b^3+3 ell|<195; and both conjugates obey
    # |a|<87/4, |b a+ell|<245/8.
    for b in list(range(-7, 0)) + list(range(1, 8)):
        for T in range(b * b - 12, b * b + 13):
            stages["b_T_box"] += 1
            root_budget_rhs = 51 - 4 * abs(T - b * b)
            if root_budget_rhs <= 0:
                continue
            # delta is positive, congruent to T^2 mod 4, nonsquare, and
            # 16 delta < (51-4|T-b^2|)^2.
            for delta in range(1, (root_budget_rhs * root_budget_rhs - 1) // 16 + 1):
                stages["b_T_delta_raw"] += 1
                if (T * T - delta) % 4 != 0 or not nonsquare(delta):
                    continue
                stages["trace_discriminant"] += 1
                N = (T * T - delta) // 4
                # Both a=(T +/- sqrt(delta))/2 have absolute value <87/4:
                # 2|T|+2sqrt(delta)<87.
                if not strict_two_embedding_bound(2 * T, 2, delta, 87):
                    continue
                stages["a_embedding"] += 1
                # Since the strict integer inequality is equivalent to
                # |b^3+3 ell|<=24, this interval is complete and minimal.
                ell_lo = -((24 + b**3) // 3)
                ell_hi = (24 - b**3) // 3
                for ell in range(ell_lo, ell_hi + 1):
                    stages["ell_raw"] += 1
                    if 8 * abs(b**3 + 3 * ell) >= 195:
                        continue
                    stages["ell_budget"] += 1
                    # Both ba+ell have |.|<245/8:
                    # 4|bT+2ell|+4|b|sqrt(delta)<245.
                    if not strict_two_embedding_bound(
                        4 * (b * T + 2 * ell), 4 * abs(b), delta, 245
                    ):
                        continue
                    stages["linked_embedding"] += 1
                    rows.append({"b": b, "T": T, "N": N, "delta": delta, "ell": ell})

    # Independent small audit of the integer radical comparison.
    assert strict_two_embedding_bound(0, 1, 2, 2)
    assert not strict_two_embedding_bound(0, 1, 4, 2)
    assert strict_two_embedding_bound(3, 2, 1, 6)
    assert not strict_two_embedding_bound(3, 2, 1, 5)

    by_b = Counter(row["b"] for row in rows)
    # The full conditional system is invariant under (b,ell)->(-b,-ell).
    # This is an inexpensive completeness guard for the integer ell interval.
    for b in range(1, 8):
        assert by_b[b] == by_b[-b], (b, by_b[b], by_b[-b])
    payload = {
        "status": "CONDITIONAL_ON_UNAUDITED_ROOT_ENVELOPE",
        "description": (
            "Exact integer enumeration conditional on the parent-supplied degree-six "
            "root/trace bounds. It is not a proof that the listed rows are all possible "
            "until those bounds are independently audited."
        ),
        "core": "P=1+b*x+a*x^2+(b*a+ell)*x^3+d*x^4+e*x^5+x^6; a^2-T*a+N=0",
        "bounds": {
            "b": "b in {-7,...,-1,1,...,7}",
            "trace_discriminant": "16*delta < (51-4*abs(T-b^2))^2, delta>0 nonsquare",
            "ell": "8*abs(b^3+3*ell)<195",
            "a_embeddings": "both |a|<87/4",
            "linked_embeddings": "both |b*a+ell|<245/8",
            "d_e_box_not_enumerated": "d in [-21,21], e in [-7,7]",
        },
        "stages": dict(stages),
        "surviving_mother_rows": len(rows),
        "estimated_d_e_expansion": len(rows) * 43 * 15,
        "by_b": {str(k): by_b[k] for k in sorted(by_b)},
        "rows": rows,
        "python": sys.version.split()[0],
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    target = HERE / "core6_linked_conditional_ledger_corrected.json"
    target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(f"status={payload['status']}")
    print("stages=" + json.dumps(payload["stages"], sort_keys=True))
    print(f"surviving_mother_rows={len(rows)}")
    print(f"estimated_d_e_expansion={payload['estimated_d_e_expansion']}")
    print("by_b=" + json.dumps(payload["by_b"], sort_keys=True))
    print(f"ledger_sha256={hashlib.sha256(target.read_bytes()).hexdigest()}")


if __name__ == "__main__":
    main()
