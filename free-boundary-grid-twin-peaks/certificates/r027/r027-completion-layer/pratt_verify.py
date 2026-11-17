#!/usr/bin/env python3
"""Pure-standard-library verifier for the recursive Lucas/Pratt certificate."""
from __future__ import annotations
from math import gcd
from pathlib import Path
import json

data = json.loads((Path(__file__).parent / "pratt-certificates.json").read_text())
nodes = data["nodes"]
done: set[int] = set()

def verify(p: int) -> None:
    if p in done: return
    node = nodes.get(str(p))
    assert node is not None and int(node["p"]) == p
    if p == 2:
        assert node == {"p": "2", "base": "two"}
        done.add(p); return
    factors = {int(q): int(e) for q, e in node["p_minus_1_factorization"].items()}
    assert p > 2 and all(q >= 2 and e >= 1 for q, e in factors.items())
    product = 1
    for q, e in factors.items():
        verify(q)
        product *= q**e
    assert product == p - 1
    a = int(node["witness"])
    assert 1 < a < p
    assert pow(a, p-1, p) == 1
    # Full factorization plus these gcd tests is Lucas's criterion: modulo
    # every prime divisor of p, a has order divisible by every q^e in p-1;
    # hence that divisor is 1 mod (p-1), so it must equal p.
    for q in factors:
        assert gcd(pow(a, (p-1)//q, p) - 1, p) == 1
    done.add(p)

for root in map(int, data["roots"]): verify(root)
print(json.dumps({"all_verified": True, "roots": data["roots"], "verified_nodes": len(done)}, indent=2))
