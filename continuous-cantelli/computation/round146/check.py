import sympy as s,json
t,h,w=s.symbols('t h w',positive=True)
shift=-t+s.sqrt(t*t+2*h*t**s.Rational(3,2))
assert s.simplify(shift+shift**2/(2*t)-h*s.sqrt(t))==0
clock=w*(1+4*h/w)**s.Rational(1,4)
print(json.dumps(dict(sympy=s.__version__,likelihood_checked=True,clock=str(s.series(clock,w,s.oo,3)))))
