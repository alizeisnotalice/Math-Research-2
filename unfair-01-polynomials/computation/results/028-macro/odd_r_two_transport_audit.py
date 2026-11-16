#!/usr/bin/env python3
"""Finite exact checks for the odd-r conditional two-transport theorem.

The scan is only an accounting cross-check in a small box.  It does not
replace the theorem's explicit injectivity and alias assumptions.
"""
from collections import Counter, defaultdict
from math import gcd
import sympy as sp


def P(r, a):
    ans = Counter({(0, i*a): 1 for i in range(r)})
    ans[(a, 0)] += 1
    ans[(a, -a)] += 1
    return ans


def mul(A, B):
    C = Counter()
    for (r,e),x in A.items():
        for (s,f),y in B.items(): C[r+s,e+f] += x*y
    return C


def aliases(Tcnt, d):
    T=set(Tcnt); im=defaultdict(list)
    for layer,e in T: im[layer*d+e].append((layer,e))
    if not all(len(v)==1 for v in im.values()): return None
    TT=set(mul(Tcnt,Tcnt))
    return [((layer,e),im[layer*d+e][0]) for layer,e in TT
            if layer*d+e in im and im[layer*d+e][0] != (layer,e)]


# r-independent row identity.
La,Lb,Ua,Ub,D,E,z=sp.symbols('La Lb Ua Ub D E z')
V=La*Lb; R=Ua*Lb+z*D; S=La*Ub-D; T=z*E; W=Ua*Ub-E
assert sp.expand(R*S*W-S**2*T-W**2*V-(La*E-Ua*D)*(Lb*W-z*Ub*S))==0

counts={}
for r in (3,5,7):
    good_pairs=admissible=0
    for a in range(2,13):
        for b in range(a+1,2*a):
            if b//gcd(a,b) < r: continue
            good_pairs += 1
            # Every component of L_a U_b has t=u^b width at most one.
            groups=defaultdict(list)
            for i in range(r):
                for j in (0,1):
                    e=i*a-j*b
                    residue=e % b
                    groups[residue].append((e-residue)//b)
            assert len(groups)==r
            assert all(max(q)-min(q)<=1 for q in groups.values())
            Tcnt=mul(P(r,a),P(r,b))
            assert set(Tcnt.values())=={1}
            # The two macro terms (a,-a) and (0,(d-1)a) collide for d<=r.
            assert all(aliases(Tcnt,d) is None for d in range(1,r+1))
            for d in range(2,41):
                data=aliases(Tcnt,d)
                if data is None: continue
                zero=[(layer,e) for layer,e in Tcnt if layer*d+e==0]
                assert zero==[(0,0)]
                if {(x[0][0],x[1][0]) for x in data} <= {(a,b),(2*a,a+b)}:
                    admissible += 1
    counts[r]=(good_pairs,admissible)

print('row syzygy exact; component-width condition verified')
print('small-box (eligible a,b pairs, macro-admissible triples):',counts)
