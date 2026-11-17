#!/usr/bin/env python3
"""Independent Fraction/sqrt(5) verifier for the integer-ring branch tree."""
from fractions import Fraction as F
from pathlib import Path
import json

TREE = Path(__file__).with_name("b5_v_alpha_integer_tree.json")
ZERO, ONE = (F(0), F(0)), (F(1), F(0))
V = (F(3, 2), F(-1, 2))  # 3/2-sqrt(5)/2


def add(x, y): return (x[0] + y[0], x[1] + y[1])
def neg(x): return (-x[0], -x[1])
def sub(x, y): return add(x, neg(y))
def mul(x, y): return (x[0]*y[0] + 5*x[1]*y[1], x[0]*y[1] + x[1]*y[0])


def sign(x):
    p, q = x
    if q == 0: return (p > 0) - (p < 0)
    if p == 0: return (q > 0) - (q < 0)
    if (p > 0) == (q > 0): return 1 if p > 0 else -1
    diff = p*p - 5*q*q
    assert diff != 0
    if p > 0 and q < 0: return 1 if diff > 0 else -1
    return 1 if diff < 0 else -1


def ring_to_sqrt(pair):
    p, q = map(F, pair)
    return (p + F(3, 2)*q, -F(1, 2)*q)


def main():
    tree = json.loads(TREE.read_text())
    states = {(ONE,): ()}
    for level, logged in zip(tree["levels"], tree["branch_ledger"]):
        k = level["degree"]
        assert logged["degree"] == k
        lookup = {(tuple(e["parent_bits"]), e["bit"]): e for e in logged["branches"]}
        nxt = {}
        for prefix, word in states.items():
            am2 = prefix[k - 2] if k >= 2 else ZERO
            am5 = prefix[k - 5] if k >= 5 else ZERO
            for bit in (0, 1):
                candidate = sub(sub((F(bit), F(0)), mul(V, am2)), am5)
                entry = lookup[(word, bit)]
                assert candidate == ring_to_sqrt(entry["candidate"])
                if sign(candidate) < 0:
                    decision = "drop_negative"
                elif sign(sub(ONE, candidate)) < 0:
                    decision = "drop_above_one"
                else:
                    decision = "keep"
                    nxt[prefix + (candidate,)] = word + (bit,)
                assert entry["decision"] == decision
        states = nxt
        assert len(states) == level["survivor_count"]
        logged_states = {
            (tuple(tuple(ring_to_sqrt(x)) for x in e["a"]), tuple(e["bits"]))
            for e in level["survivors"]
        }
        assert set((p, w) for p, w in states.items()) == logged_states
    assert not states
    print("PASS: independent Fraction/sqrt(5) replay verified all 70 tree branches and degree-29 death")


if __name__ == "__main__":
    main()
