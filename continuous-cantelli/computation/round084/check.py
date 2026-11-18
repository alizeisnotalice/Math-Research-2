"""Symbolic endpoint quadratic and quartic backward heat test."""
import sympy as s
x,A,B,u,v=s.symbols('x A B u v')
p=(x+A)*(B-x)
def heat(poly,t):
    return s.expand(sum((t/2)**j*s.diff(poly,x,2*j)/s.factorial(j) for j in range(s.degree(poly,x)//2+1)))
for poly in [p, p*(x-v)**2, p**2, (x-v)**6]:
    assert s.expand(heat(heat(poly,-u),u)-poly)==0
print(heat(-p,-u))
print('4 exact inverse heat identities')
