import json
from math import comb, factorial
import sympy as s

a,z,x = s.symbols('a z x', positive=True)
p = s.Poly(1,a)
rows = []
for n in range(13):
    co = [int(p.nth(i)) for i in range(n+1)]
    assert co[-1] == 1 and min(co)>0
    kernel = s.exp(-z*x)*p.as_expr().subs(a,z*x)/z**(2*n+1)
    nxt = s.Poly((2*n+1+a)*p.as_expr()-a*p.diff().as_expr(),a)
    target = s.exp(-z*x)*nxt.as_expr().subs(a,z*x)/z**(2*n+3)
    assert s.expand(((-s.diff(kernel,z)/z-target)*s.exp(z*x))).simplify()==0
    rows.append(dict(n=n, integral=sum(c*factorial(i) for i,c in enumerate(co)),coefficients=co))
    p=nxt
print(json.dumps(dict(sympy=s.__version__,rows=rows),indent=2))
