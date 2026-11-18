"""Exact SymPy coefficient verification after blowups at the heat and space endpoints."""
import json
from pathlib import Path
import sympy as s

p,b,w,v=s.symbols('p beta w v', real=True)
scaled=s.exp(2*b*w*p-w*w*p**2)/(1-b*p+w*p*p)**2*(1+2*v*p*p/(1-b*p+w*p*p)**2)**(-s.Rational(3,2))
polynomial=s.Poly(s.series(scaled,p,0,4).removeO().expand(),w,v)
series=s.expand(sum(c*s.factorial(i)*s.factorial(j) for (i,j),c in polynomial.terms()))
assert series.coeff(p,1)==4*b
assert series.coeff(p,2)==11*b*b-7
assert series.subs(b,0)==1-7*p*p
C,L=s.symbols('C L')
R=b/p+L-b*b-C
D=s.series(2*p*p*R/(1-2*p*p*R),p,0,3).removeO().expand()
assert D.coeff(p,2)==2*L+2*b*b-2*C
Path('symbolic.json').write_text(json.dumps(dict(normalized_endpoint_series=str(series),inverse_through_second=str(D),sympy=s.__version__),indent=2)+'\n')
