"""Exact macro-coordinate audit for H_d(x^2)H_d(x^3).

All coordinates are (v-layer, u-offset), representing u**offset*v**layer
after v=u**d.  This checks only finite support arithmetic; the accompanying
proof supplies the general nonnegative-factor arguments.
"""
import sympy as sp
import json
from collections import Counter


def target_coords():
    return (
        {(0, e) for e in (0, 2, 3, 4, 5, 6, 7, 8, 10)}
        | {(2, e) for e in (-2, 0, 1, 3, 4, 6)}
        | {(3, e) for e in (-3, -1, 0, 1, 2, 4)}
        | {(5, e) for e in (-5, -3, -2, 0)}
    )


def actual_aliases(d):
    target = target_coords()
    image = {e + d * r: (r, e) for r, e in target}
    assert len(image) == len(target) == 25
    actual = {(r + s, e + f) for r, e in target for s, f in target}
    return sorted(
        ((r, e), image[e + d * r])
        for r, e in actual
        if e + d * r in image and (r, e) != image[e + d * r]
    )


def template(a):return [(0,0),(0,a),(0,2*a),(a,-a),(a,0)]
macro=Counter((r+s,e+f) for r,e in template(2) for s,f in template(3))
assert set(macro)==target_coords() and set(macro.values())=={1}

for d in range(13, 20):
    assert all(
        (actual[0], canonical[0]) in {(2, 3), (4, 5)}
        for actual, canonical in actual_aliases(d)
    ), (d, actual_aliases(d))
for d in (18, 20, 21, 22, 23, 24, 25):
    assert not actual_aliases(d), (d, actual_aliases(d))

L2, U2, L3, U3, D, E, z = sp.symbols("L2 U2 L3 U3 D E z")
V = L2 * L3
R = U2 * L3 + z * D
S = L2 * U3 - D
T = z * E
W = U2 * U3 - E
assert sp.expand(
    R * S * W - S**2 * T - W**2 * V
    - (L2 * E - U2 * D) * (W * L3 - z * S * U3)
) == 0

expected={11: [[[0, 20], [2, -2]], [[2, 8], [3, -3]], [[2, 10], [3, -1]], [[2, 11], [3, 0]], [[2, 12], [3, 1]], [[2, 13], [3, 2]], [[4, 6], [5, -5]], [[4, 8], [5, -3]], [[4, 9], [5, -2]]], 12: [[[2, 9], [3, -3]], [[2, 11], [3, -1]], [[2, 12], [3, 0]], [[2, 13], [3, 1]], [[2, 14], [3, 2]], [[2, 16], [3, 4]], [[4, 7], [5, -5]], [[4, 9], [5, -3]], [[4, 10], [5, -2]], [[4, 12], [5, 0]]]}
for d in (11,12):
 assert [[list(x),list(y)] for x,y in actual_aliases(d)]==expected[d]
for d in [11,*range(13,20)]:
 assert all(canonical!=(5,0) for actual,canonical in actual_aliases(d))
u=sp.symbols('u');L2u=1+u**2+u**4;L3u=1+u**3+u**6
assert sp.gcd(L2u,1+u**2)==1 and sp.gcd(L3u,1+u**3)==1
assert sp.expand(L3u.subs(u,u+1))==u**6+6*u**5+15*u**4+21*u**3+18*u**2+9*u+3
for d in range(18):
 H=lambda k:1+u**k+u**(2*k)+u**(k*(d+17))+u**(k*(d+18))
 assert sp.gcd(L2u*L3u,H(2)*H(3))==1
rows=[{'d':d,'aliases':[[list(x),list(y)] for x,y in actual_aliases(d)]} for d in range(11,26)]
print(json.dumps({'status':'passed_exact_macro_rows_and_identity','rows':rows,'row_identity_verified':True,'gcd_crosscheck_residue_classes':18,'scope':'Finite coordinate and algebra audits; infinite d uses the accompanying proof.'},indent=2))
