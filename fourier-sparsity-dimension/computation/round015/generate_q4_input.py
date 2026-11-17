#!/usr/bin/env python3
"""Exact GF(2) construction of the three Q=4 selector-plus-Cdiag codes."""
import hashlib, json
from pathlib import Path

N = 64

def rref(vectors):
    pivots = {}
    for v in vectors:
        while v:
            p = v.bit_length() - 1
            if p in pivots: v ^= pivots[p]
            else: pivots[p] = v; break
    return [pivots[p] for p in sorted(pivots, reverse=True)]

def mask(pred): return sum(pred(q) << q for q in range(N))
def bit(q, j): return (q >> j) & 1

def make(a_name):
    def a(q):
        x, y = bit(q, 0), bit(q, 2)
        return {"x": x, "y": y, "x_plus_y": x ^ y}[a_name]
    gens = []
    # h_{w0,w1} in span{1,a}, independently on each of four w cells.
    for w0 in range(2):
        for w1 in range(2):
            for e in (lambda q: 1, a):
                gens.append(mask(lambda q, w0=w0, w1=w1, e=e: int(bit(q,3)==w0 and bit(q,5)==w1) * e(q)))
    # z0*E.
    for e in (lambda q: 1, a):
        gens.append(mask(lambda q, e=e: bit(q,1) * e(q)))
    # z1*a_{w0} in E, independently for w0=0,1.
    for w0 in range(2):
        for e in (lambda q: 1, a):
            gens.append(mask(lambda q, w0=w0, e=e: bit(q,4) * int(bit(q,3)==w0) * e(q)))
    s_basis = rref(gens)
    if len(s_basis) != 14: raise AssertionError((a_name, len(s_basis)))
    cdiag = [mask(lambda q: 1), mask(lambda q: bit(q,0)), mask(lambda q: bit(q,2))]
    basis = rref(s_basis + cdiag)
    if len(basis) != 15: raise AssertionError((a_name, len(basis)))
    return {"name": a_name, "rank": len(basis), "basis_masks": basis, "E_rows": [],
            "S_rank": len(s_basis), "Cdiag_rank": 3, "E_basis": ["1", a_name]}

out = {"schema": "r015-q4-selector-affine-v1", "domain_bits": "x,z0,y,w0,z1,w1",
       "base_phase": "x*y + z0*w0 + z1*w1", "records": [make(a) for a in ("x", "y", "x_plus_y")]}
text = json.dumps(out, indent=2, sort_keys=True) + "\n"
Path(__import__('sys').argv[1]).write_text(text)
print(json.dumps({"records": len(out["records"]), "ranks": [r["rank"] for r in out["records"]],
                  "sha256": hashlib.sha256(text.encode()).hexdigest()}))
