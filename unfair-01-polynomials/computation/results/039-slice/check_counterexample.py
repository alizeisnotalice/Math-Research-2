#!/usr/bin/env python3
"""Exact algebraic audit of a non-rational section fixed point.

p=2,q=3; a is the positive root a^3-a^2-1=0.  The example is deliberately
only a boundary example: it satisfies binary H,U and polynomial F,G, but its
full crossed products/C are audited separately and fail.
"""
from pathlib import Path
import hashlib, json, sys
import sympy as S

x,y,a=S.symbols('x y a')
rel=a**3-a**2-1
def red(z):
    return S.rem(S.Poly(S.expand(z),a),S.Poly(rel,a)).as_expr()
def coeffs(P):
    p=S.Poly(S.expand(P),x)
    return {i:red(v) for (i,),v in p.terms()[::-1] if red(v)!=0}
def binary(P):
    return set(coeffs(P).values())=={S.Integer(1)}

R=1+a*y
U=1+x**2+x**6
F=1+(1-a)*x**2+(a**2-a)*x**4
G=(1+x)*(1+a*x**2)
H=S.expand(U*(1+x))
A=S.expand(F.subs(x,x**2)*G.subs(x,x**3))
B=S.expand(F.subs(x,x**3)*G.subs(x,x**2))
C=S.expand(H.subs(x,x**2)*H.subs(x,x**3))
assert red(S.expand(F*(1+a*x**2)-U))==0
assert red(S.expand(F*G-H))==0
assert binary(U) and binary(H)
# exact section of G at q=3, with y reindexing
RG=sum(red(v)*y**(i//3) for i,v in coeffs(G).items() if i%3==0)
assert red(RG.subs(y,y)-R)==0
poly=S.Poly(rel,a)
isolating=S.intervals(poly,eps=S.Rational(1,10)**12)
positive=[z for z in isolating if z[0][0]>0]
assert len(positive)==1
apos=positive[0][0]
# a>1 makes the x^4 coefficient of A exactly 1-a negative.
assert coeffs(A)[4]==1-a
record={
 'status':'passed_exact_boundary_counterexample',
 'parameters':{'p':2,'q':3,'minimal_polynomial':'a^3-a^2-1','positive_isolating_interval':[str(apos[0]),str(apos[1])]},
 'polynomials':{'R(y)':str(R),'U':str(U),'F':str(F),'G':str(G),'H':str(H)},
 'identities':{'U=F R(x^p)':True,'H=F G':True,'R=section_q(G)':True,'T=H/U':'1+x'},
 'integrality':{'a_rational':False,'R_not_in_Qy':True,'F_not_in_Qx':True,'G_not_in_Qx':True},
 'full_cross_failure':{
  'A_x4':str(coeffs(A)[4]),
  'reason':'a>1, hence 1-a<0; A is not coefficientwise nonnegative',
  'C_binary':binary(C),
  'C_coefficients_not_01':{str(i):str(v) for i,v in coeffs(C).items() if v not in (0,1)}
 },
 'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'python':sys.version,'sympy':S.__version__,
 'scope':'Exact counterexample to fixed-point uniqueness/rationality from H,U binary plus polynomial F,G alone; it is not a complete crossed-model or original-conjecture counterexample.'
}
Path(__file__).with_suffix('.json').write_text(json.dumps(record,indent=2,sort_keys=True)+'\n')
print(json.dumps(record,indent=2,sort_keys=True))
