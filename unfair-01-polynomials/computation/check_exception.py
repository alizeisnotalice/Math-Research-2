#!/usr/bin/env python3
"""Check the six complete diagonals and exact positive residual certificate."""
import json
import sympy as S

x,p,q,r,s,v,w,u = S.symbols('x p q r s v w u')
A=1+p*x+r*x**2+s*x**4+x**9
B=1+q*x+u*x**3+v*x**4+w*x**6+x**11
C=S.Poly(A*B,x)
expected={1:p+q,2:r+p*q,6:w+r*v,10:s*w+q,13:r+v,15:s+w}
assert all(S.expand(C.nth(k)-value)==0 for k,value in expected.items())
e={k:S.expand(value-1) for k,value in expected.items()}
rhs=e[1]-e[10]+w*(e[15]-e[6]+r*(e[13]-e[2]))
target=p*(1-r*q*w)
assert S.expand(rhs-target)==0
assert S.expand((target-p**2).subs(q,1-p)-p*(1-p)*(1-r*w))==0
print(json.dumps({'I':[0,1,2,4,9],'J':[0,1,3,4,6,11],
                  'diagonals':{k:str(z) for k,z in expected.items()},
                  'identity_verified':True,
                  'positive_expression':'p*(1-r*q*w) >= p^2 > 0',
                  'assumptions':['p>0','q>0','p+q=1','0<r<=1','0<w<=1'],
                  'sympy':S.__version__},indent=2))
