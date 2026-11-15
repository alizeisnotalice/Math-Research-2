#!/usr/bin/env python3
import json
import sympy as S
s,t=S.symbols('s t',real=True)
checks=[]
for r in range(2,13):
    beta=[S.Integer(1)]
    for k in range(1,r+2):
        beta.append(S.expand(1-s*beta[k-1]-(t*beta[k-r] if k>=r else 0)))
    target=(1-(-s)**(r+2))/(1+s)+(2*s-1)*t
    assert S.cancel(beta[r+1]-target)==0
    assert S.cancel((target-(2*s-1)*t+2*s-1)-(s+2*s**2-(-s)**(r+2))/(1+s))==0
    x=S.Symbol('x')
    for d in (1,2,3):
        aa=1+s*x**d+t*x**(r*d)+x**((2*r+1)*d)
        bs=S.symbols('b0:'+str((2*r+1)*d+1))
        bb=sum(b*x**k for k,b in enumerate(bs))
        poly=S.Poly(aa*bb,x)
        assert S.expand(poly.nth((r+1)*d)-(bs[(r+1)*d]+s*bs[r*d]+t*bs[d]))==0
        assert S.expand(poly.nth((2*r+1)*d)-(bs[(2*r+1)*d]+s*bs[2*r*d]+t*bs[(r+1)*d]+bs[0]))==0
    checks.append({'r':r,'forced_beta_rplus1':str(beta[r+1]),'identity_verified':True,'d_checked':[1,2,3]})
print(json.dumps({'sympy':S.__version__,'cases':checks,'general_proof':'proofs/005-unbounded-four-term-family.md'},indent=2))
