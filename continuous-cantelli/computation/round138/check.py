from sympy import symbols, diff, exp, simplify, sqrt, Rational
import sympy,json
s,d,p,w=symbols('s d p w',positive=True)
A=d**2/2
B=s**(-p-Rational(1,2))*exp(A/s)
z=d/s+w/sqrt(s)
assert simplify(z*d-s*z**2/2-A/s+w**2/2)==0
assert simplify(-diff(B,s)/B-A/s**2-(p+Rational(1,2))/s)==0
print(json.dumps(dict(sympy=sympy.__version__,completed_square=True,layer_derivative=True)))
