#!/usr/bin/env python3
import json
import sympy as S
x,s,t,u,v,w=S.symbols('x s t u v w')
expected={3:s+u,5:t+v,10:u+t*v,12:v+s*w}
for extra in (False,True):
    I=[0,3,5,7];J=[0,3,5,6,8,9,11,15,17]+([13] if extra else [])
    a={0:1,3:s,5:t,7:1}
    b={j:S.Symbol('b'+str(j)) for j in J};b.update({0:1,17:1,3:u,5:v,9:w})
    C=S.Poly(sum(a[i]*x**i for i in I)*sum(b[j]*x**j for j in J),x)
    assert all(S.expand(C.nth(k)-z)==0 for k,z in expected.items())
e={k:z-1 for k,z in expected.items()}
assert S.expand(e[5]-e[12]+w*(e[3]-e[10])-t*(1-v*w))==0
assert S.expand((t*(1-v*w)-t**2).subs(v,1-t)-t*(1-t)*(1-w))==0
print(json.dumps({'sympy':S.__version__,'diagonals':{k:str(z) for k,z in expected.items()},
                  'extra_b13_checked':True,'identity':'e5-e12+w*(e3-e10) = t*(1-v*w) >= t^2 > 0',
                  'identity_verified':True},indent=2))
