#!/usr/bin/env python3
"""Exact symbolic identities for the all-even-k>=8 linear threshold."""
import json
from fractions import Fraction as F
import sympy as S

s,k,y,z,u,v,H=S.symbols('s k y z u v H')
h=1+s;a=2*s-1;L=k+2+(k+1)*s;w=(1+L*y)/h**2
w1=w-(k+2)*y;w2=w1+(k+3)*s*y
q=(2*s*s-s)*h*h*w*w+(1-z)*(1-s-2*s*s*z)*w+(1-z)**2*(s*w1-w2)
qu=s*a*H*H-a*H*u+(a-3*s*s*y)/h*u*u
assert S.cancel(q-qu.subs({H:h*w,u:1-z}))==0
disc=S.discriminant(qu,u)
assert S.cancel(disc-a*H*H/h*(-a*(3*s-1)+12*s**3*y))==0
f=a*(3*s-1)/(12*s**3)
assert S.cancel(S.diff(f,s)+(6*s*s-10*s+3)/(12*s**4))==0
assert F(3,4)**9<F(3,32)
L8=10+9*s;c=3*s*s/L8
lower=s*a*H*H-(a*u+c*u*u)*H+(a+c)/h*u*u
assert S.cancel((qu-lower).subs(H,(1+L*y)/h)-3*s*s*y*u*u/h*(L/L8-1))==0
assert S.cancel(S.discriminant(lower,H)-u*u*((a+c*u)**2-4*s*a*(a+c)/h))==0
T=4*s*a*L8*(a*L8+3*s*s)-h*(a*L8+6*s*s)**2
expected=936*v**5+4122*v**4+5856*v**3+2912*v*v+S.Rational(15395,32)*v+S.Rational(1435,128)
assert S.expand(T.subs(s,v+S.Rational(3,4))-expected)==0
E8=F(2,3)**8*F(48,3)+F(5,3)*F(2,3)**19
assert E8==F(2179403776,3486784401)<1
assert S.expand(2*(5*k+13)-3*(5*k+8))==2-5*k
print(json.dumps({'verified':True,'sympy':S.__version__,'small_parameter_bound':str(E8),
                  'scope':'all even k>=8 and r>=k+3; symbolic identities plus analytic sign proofs'}))
