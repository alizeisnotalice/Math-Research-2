#!/usr/bin/env python3
"""Rebuild the nine-point support equations and check their exact identity."""

from collections import defaultdict
from hashlib import sha256
import json

import sympy as sp


I = (0,2,3,4,7,8,9,14,20)
J = (0,2,3,5,8,9,10,15,21)
representations = defaultdict(list)
for i in I:
    for j in J:
        representations[i+j].append((i,j))

needed = {
    0: [(0,0)],
    2: [(0,2),(2,0)],
    3: [(0,3),(3,0)],
    4: [(2,2),(4,0)],
    5: [(0,5),(2,3),(3,2)],
    20: [(20,0)],
    21: [(0,21)],
    25: [(4,21),(20,5)],
}
assert all(representations[k] == expected for k,expected in needed.items())
aa = {i:sp.Symbol(f"a{i}") for i in I}
bb = {j:sp.Symbol(f"b{j}") for j in J}
normalize = {aa[0]:1,bb[0]:1,aa[20]:1,bb[21]:1}
residuals = {k:sum(aa[i]*bb[j] for i,j in representations[k]).subs(normalize)-1
             for k in (2,3,4,5,25)}
x,u,y,v = aa[2],bb[2],aa[3],bb[3]
positive = u**2*v+x**2*y
certificate = ((1-x*u)*residuals[3]
               + ((u+1)*v+(x+1)*y)*residuals[2]
               - residuals[5]+residuals[25]-residuals[4])
assert sp.expand(positive-certificate) == 0

out = dict(status="verified", supports=[I,J], required_rows=needed, convolution_rows=len(representations),
           identity=str(positive), expanded_difference=str(sp.expand(positive-certificate)),
           sympy=sp.__version__)
out["input_sha256"] = sha256(json.dumps([I,J],separators=(",",":")).encode()+b"\n").hexdigest()
print(json.dumps(out,sort_keys=True))
