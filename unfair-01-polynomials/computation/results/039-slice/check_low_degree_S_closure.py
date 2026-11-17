#!/usr/bin/env python3
"""Exact regression checks for the low-degree S closure lemma.

The mathematical all-b proof is in low_degree_S_closure.md.  This script
checks its symbolic branch inequalities with rationals and independently
tests the residue-class divisibility identity on exhaustive small examples.
"""
from pathlib import Path
import hashlib,json,sys
from itertools import product
import sympy as S
x,y=S.symbols('x y')

def sections(poly,p):
    P=S.Poly(poly,x);out=[]
    for r in range(p):
        out.append(sum(c*y**((k-r)//p) for (k,),c in P.terms() if k%p==r))
    return out

# Exhaustive small polynomial identity: H=Q(x^p)W implies Q|all p-sections.
divisibility=[]
for p in (2,3,4):
    Q=1+y+y*y
    for bits in product((0,1),repeat=5):
        W=sum(bits[i]*x**i for i in range(5))
        H=S.expand(Q.subs(y,x**p)*W)
        for R in sections(H,p):
            assert S.rem(S.Poly(R,y),S.Poly(Q,y))==0
        divisibility.append((p,bits))

# Exact parameter branch audit for Q=y^2+b y+epsilon.
# epsilon=-1 always has opposite-sign roots.  epsilon=1, disc>0 means
# integer b>=3 or b<=-3; the b>=3 branch has Q(-1/2)<0<Q(0).
branch={
 'constant_minus_one': {'product_of_roots':-1,'conclusion':'opposite signs, hence a positive root'},
 'constant_one_positive_roots': {'condition':'b<=-3','sum_of_roots':'-b>0','product':1},
 'constant_one_negative_small_root': {'condition':'b>=3','Q(-1/2)':'5/4-b/2<0','Q(0)':1,
                                      'conclusion':'a root lies strictly in (-1/2,0)'}
}
for b in range(3,1001):
    assert S.Rational(5,4)-S.Rational(b,2)<0
for b in range(-1000,-2):
    # both roots positive: sum -b and product 1, discriminant positive.
    assert -b>0 and b*b-4>0

# finite geometric bound is exact for every finite degree d at radius one-half.
for d in range(1,200):
    assert sum(S.Rational(1,2)**j for j in range(1,d+1))<1
record={'status':'passed_exact_low_degree_closure_regressions','divisibility_examples':len(divisibility),
 'quadratic_parameter_branches':branch,'finite_geometric_bound_checked_through_degree':199,
 'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'python':sys.version,'sympy':S.__version__,
 'scope':'Exact regression of the symbolic residue and quadratic branches; the unbounded theorem is supplied by the analytic proof, not a finite loop.'}
Path(__file__).with_suffix('.json').write_text(json.dumps(record,indent=2,sort_keys=True)+'\n');print(json.dumps(record,indent=2,sort_keys=True))
