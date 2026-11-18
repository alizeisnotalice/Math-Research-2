#!/usr/bin/env python3
"""Directed interval certificates for one coarse failure and one fine coupling."""
import json
from fractions import Fraction
from pathlib import Path

import mpmath as mp

iv=mp.iv
iv.dps=45
F=Fraction
depth=6
subboxes=128
sigma=iv.mpf(7)/2


def e(x):
    return iv.mpf(x.numerator)/iv.mpf(x.denominator)


def box(a,b):
    return iv.mpf([e(a).a,e(b).b])


def prefix(n):
    p=[0]
    for _ in range(n):
        p=[k for old in p for k in (3*old,3*old+2)]
    return [box(F(i,3**n),F(i+1,3**n)) for i in p]


def kernel(a,x):
    z=a-x
    return iv.exp(-((z/sigma)**2+a*(a+4)+x*(x+4))/2)*iv.cos(2*iv.pi*z)


leaves=prefix(depth)
gap=[box(F(subboxes+i,3*subboxes),F(subboxes+i+1,3*subboxes)) for i in range(subboxes)]
E=sum((kernel(a,x) for a in leaves for x in leaves),iv.mpf(0))/len(leaves)**2
G=sum((kernel(a,x) for a in leaves for x in gap),iv.mpf(0))/(len(leaves)*len(gap))
assert E.a>iv.mpf(".060")
assert G.b<iv.mpf("-.045")
L=(2+iv.exp(-iv.mpf('0.5'))/sigma+2*iv.pi)/(3*3**4)
assert L.b<iv.mpf(".035")
assert (E-L).a>iv.mpf(".025")
result=dict(m=0,b="49/4",D="49*pi/2",prefix_depth=depth,uniform_gap_boxes=subboxes,
            energy_normalized=str(E),true_gap_j0_real_normalized=str(G),
            j4_lipschitz_bound=str(L),j4_real_strict_lower_bound=str(E-L))
Path("certificate.json").write_text(json.dumps(result,indent=2)+"\n")
for k,v in result.items():
    print(k,v)
