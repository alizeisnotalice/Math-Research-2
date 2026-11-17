#!/usr/bin/env python3
"""Exact all-prefix tree for B=1+(2-sqrt(3))*x^2+x^5.

Coefficients are Fraction pairs p+q*sqrt(3), not the parent's Z[v] floor
representation.  Every dropped branch records the exact sign witness.
"""
from fractions import Fraction as F
from pathlib import Path
import json

ZERO, ONE = (F(0), F(0)), (F(1), F(0))
V = (F(2), F(-1))
OUT = Path(__file__).with_name("b5_v_2_minus_sqrt3_fraction_tree.json")


def add(x, y): return (x[0] + y[0], x[1] + y[1])
def neg(x): return (-x[0], -x[1])
def sub(x, y): return add(x, neg(y))
def mul(x, y): return (x[0]*y[0] + 3*x[1]*y[1], x[0]*y[1] + x[1]*y[0])


def sign(x):
    """Exact sign of p+q sqrt(3), using only rational square comparisons."""
    p, q = x
    if q == 0: return (p > 0) - (p < 0)
    if p == 0: return (q > 0) - (q < 0)
    if (p > 0) == (q > 0): return 1 if p > 0 else -1
    diff = p*p - 3*q*q
    assert diff != 0
    if p > 0 and q < 0: return 1 if diff > 0 else -1
    return 1 if diff < 0 else -1


def jp(x): return [str(x[0]), str(x[1])]


def main():
    states = {(ONE,): ()}
    levels, ledger = [], []
    for k in range(1, 50):
        nxt, rows = {}, []
        for prefix, word in states.items():
            am2 = prefix[k - 2] if k >= 2 else ZERO
            am5 = prefix[k - 5] if k >= 5 else ZERO
            for bit in (0, 1):
                candidate = sub(sub((F(bit), F(0)), mul(V, am2)), am5)
                lo, hi = sign(candidate), sign(sub(ONE, candidate))
                row = {"parent_bits": list(word), "bit": bit, "candidate": jp(candidate)}
                if lo < 0:
                    row.update(decision="drop_negative", witness=jp(candidate))
                elif hi < 0:
                    row.update(decision="drop_above_one", witness=jp(sub(ONE, candidate)))
                else:
                    child = prefix + (candidate,)
                    nxt[child] = word + (bit,)
                    row.update(decision="keep", child_bits=list(word + (bit,)))
                rows.append(row)
        states = nxt
        levels.append({
            "degree": k, "survivor_count": len(states),
            "survivors": [
                {"bits": list(word), "a": [jp(x) for x in prefix]}
                for prefix, word in states.items()
            ],
        })
        ledger.append({"degree": k, "branches": rows})
    assert levels[-1]["survivor_count"] == 0
    OUT.write_text(json.dumps({
        "field": "Q(sqrt(3)); v=2-sqrt(3)",
        "recurrence": "a_k=c_k-v*a_{k-2}-a_{k-5}; c_k in {0,1}; a_negative=0",
        "root": {"a": [jp(ONE)], "bits": []},
        "levels": levels, "branch_ledger": ledger,
        "claim": "No nonnegative [0,1] prefix reaches degree 49; every finite A would extend by zero and therefore cannot exist.",
    }, indent=2) + "\n")
    print("PASS: Fraction/sqrt(3) tree counts:", [x["survivor_count"] for x in levels])
    print("PASS: degree 49 has no surviving prefix; artifact:", OUT)


if __name__ == "__main__":
    main()
