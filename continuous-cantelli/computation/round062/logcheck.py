import sympy as s
import json
t,b,h=s.symbols('t b h',positive=True)
a=s.symbols('a',real=True)
P=s.series((1+t)**a,t,0,5).removeO()
assert all(s.simplify(P.coeff(t,j)-s.prod(a-k for k in range(j))/s.factorial(j))==0 for j in range(5))
# Integral of x^b log(1/x)^k is obtained by differentiating the power integral.
for k in range(5):
    prim=(-1)**k*s.diff(h**(b+1)/(b+1),b,k)
    assert s.simplify(s.diff(prim,h)-h**b*(-s.log(h))**k)==0
print(json.dumps(dict(sympy=s.__version__,real_power_orders=list(range(5)),log_degrees=list(range(5)))))
