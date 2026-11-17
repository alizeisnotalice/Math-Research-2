#!/usr/bin/env python3
"""Exact exhaustive prefix tree in Z[v], v^2=3v-1, v=(3-sqrt(5))/2.

This is intentionally independent of the earlier Fraction/sqrt(5) search.
Every coefficient is an integer pair [p,q] for p+q*v.  A rejected branch
stores the exact linear form witnessing either a_k<0 or 1-a_k<0.
"""
from pathlib import Path
import json

ZERO, ONE = (0, 0), (1, 0)
OUT = Path(__file__).with_name("b5_v_alpha_integer_tree.json")


def add(x, y): return (x[0] + y[0], x[1] + y[1])
def neg(x): return (-x[0], -x[1])
def sub(x, y): return add(x, neg(y))


def vmul(x):
    """v*(p+qv)=-q+(p+3q)v, using v^2=3v-1."""
    p, q = x
    return (-q, p + 3 * q)


def sign(x):
    """Exact sign of p+qv, v=(3-sqrt(5))/2, by integer square comparison."""
    p, q = x
    # 2*(p+qv) = A + B sqrt(5).
    A, B = 2 * p + 3 * q, -q
    if B == 0:
        return (A > 0) - (A < 0)
    if A == 0:
        return (B > 0) - (B < 0)
    if (A > 0) == (B > 0):
        return 1 if A > 0 else -1
    diff = A * A - 5 * B * B
    assert diff != 0  # Otherwise sqrt(5) would be rational, except x=0.
    if A > 0 and B < 0:
        return 1 if diff > 0 else -1
    return 1 if diff < 0 else -1


def json_pair(x): return [x[0], x[1]]


def main():
    # Map a full surviving coefficient prefix to its forced bit word c_1,...,c_k.
    states = {(ONE,): ()}
    levels, ledger = [], []
    for k in range(1, 30):
        nxt = {}
        entries = []
        for prefix, word in states.items():
            am2 = prefix[k - 2] if k >= 2 else ZERO
            am5 = prefix[k - 5] if k >= 5 else ZERO
            for bit in (0, 1):
                candidate = sub(sub((bit, 0), vmul(am2)), am5)
                s_lo, s_hi = sign(candidate), sign(sub(ONE, candidate))
                entry = {
                    "parent_bits": list(word), "bit": bit,
                    "candidate": json_pair(candidate),
                }
                if s_lo < 0:
                    entry.update(decision="drop_negative", witness=json_pair(candidate))
                elif s_hi < 0:
                    entry.update(decision="drop_above_one", witness=json_pair(sub(ONE, candidate)))
                else:
                    new = prefix + (candidate,)
                    nxt[new] = word + (bit,)
                    entry.update(decision="keep", child_bits=list(word + (bit,)))
                entries.append(entry)
        states = nxt
        levels.append({
            "degree": k,
            "survivor_count": len(states),
            "survivors": [
                {"bits": list(word), "a": [json_pair(x) for x in prefix]}
                for prefix, word in states.items()
            ],
        })
        ledger.append({"degree": k, "branches": entries})
    assert levels[-1]["survivor_count"] == 0
    payload = {
        "ring": "Z[v]/(v^2-3v+1), v=(3-sqrt(5))/2",
        "recurrence": "a_k=c_k-v*a_{k-2}-a_{k-5}; c_k in {0,1}; a_negative=0",
        "root": {"a": [[1, 0]], "bits": []},
        "levels": levels,
        "branch_ledger": ledger,
        "claim": "No nonnegative formal coefficient prefix reaches degree 29; hence no finite normalized A can satisfy A*(1+v*x^2+x^5) binary.",
    }
    OUT.write_text(json.dumps(payload, indent=2) + "\n")
    print("PASS: integer-ring tree counts:", [x["survivor_count"] for x in levels])
    print("PASS: degree 29 has no surviving prefix; artifact:", OUT)


if __name__ == "__main__":
    main()
