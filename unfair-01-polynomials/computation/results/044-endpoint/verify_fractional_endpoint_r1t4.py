#!/usr/bin/env python3
"""Exact symbolic checks for the (r,t)=(1,4) endpoint obstruction."""
from __future__ import annotations

import json
import platform
import sys

import sympy as sp

u = sp.symbols("u", real=True)
v = 1 - u
s = 1 - u*v
L = u*(u-v) + s
gap = sp.factor(1 - u + u**2*s - L**2 / 4)
expected = (1-u)*(5*u**3 - 3*u**2 + 3*u + 3) / 4
assert sp.expand(gap - expected) == 0
assert sp.expand(5*u**3 - 3*u**2 + 3*u + 3 - (5*u**3 + 3*u*(1-u) + 3)) == 0

def pairs(d):
    return [(r, d-s) for r in range(1, d) for s in range(r, d)]

actual = pairs(5)
remaining = [(2, 3), (3, 2)]
forbidden = []
for r, t in actual:
    bad = r == t or (r, t) in {(1, 2), (2, 1), (1, 3), (3, 1), (1, 4), (4, 1)}
    if not bad:
        forbidden.append((r, t))
assert forbidden == remaining

print(json.dumps({
    "python": sys.version.split()[0],
    "platform": platform.platform(),
    "sympy": sp.__version__,
    "gap_factor": str(gap),
    "degree_5_pairs": actual,
    "surviving_pairs_after_endpoint_obstructions": forbidden,
    "status": "all exact assertions passed",
}, ensure_ascii=False, indent=2))
