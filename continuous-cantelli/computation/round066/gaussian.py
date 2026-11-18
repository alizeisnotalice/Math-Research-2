import sympy as s
import json
x=s.symbols('x',real=True)
checks=0
for k in range(1,10):
    for j in range(k):
        for l in range(4):
            n=k+2*l
            poly=s.Poly(x**j*s.hermite_prob(n,x),x)
            value=sum(c*s.factorial2(m-1) for (m,),c in poly.terms() if m%2==0)
            assert value==0
            checks+=1
print(json.dumps(dict(sympy=s.__version__,gaussian_moments_checked=checks)))
