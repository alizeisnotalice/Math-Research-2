#!/usr/bin/env python3
"""Exact threshold computation and symbolic finite-error identity checks."""
from fractions import Fraction as Q
import json
import sympy as S


def threshold(k):
    assert k>=4 and k%2==0
    M=(k+3)*(k+4)//2
    delta=Q(1,4*k*(k+1))
    L=Q(2,3)**(k+1)*Q(5*k-12,36)
    p=1
    while 2**p<=11*M/L:p+=1
    R0=8*M*delta.denominator**(k+1)+10
    R=max(R0,2*(k+3)*p)
    assert Q(2*M,R0-3)<delta**(k+1)/4
    assert delta**(k+1)/4<delta**(k-1)*Q(k-2,3)
    assert Q(11*M,2**p)<L
    return {'k':k,'M':M,'delta':str(delta),'L':str(L),'p':p,'R0':R0,'R':R}


s,z,w,v = S.symbols('s z w v')
h=1+s
q=(2*s*s-s)*h*h*w*w+(1-z)*(1-s-2*s*s*z)*w+(1-z)**2*v
C1=-(1-s+2*s*s)*w-2*v
C2=2*s*s*w+v
assert S.expand(q-q.subs(z,0)-z*C1-z*z*C2)==0
# Reconstruct the fixed contribution directly from finite geometric sums.
w1,w2=S.symbols('w1 w2')
t=(1-z)/(h*w)
b1=(1+s*z)/h-t*w1
b2=(1-s*s*z)/h-t*w2
P=s*(s-t*b1)+t*b2+1-s+s*s
assert S.cancel((P-1)*h*h*w*w-q.subs(v,s*w1-w2))==0
print(json.dumps({'symbolic_identities_verified':2,'sympy':S.__version__,
                  'thresholds':[threshold(k) for k in (4,6,8,10,20)],
                  'scope':'Examples of a formula proved for every even k>=4; no exhaustive search up to R.'},indent=2))
