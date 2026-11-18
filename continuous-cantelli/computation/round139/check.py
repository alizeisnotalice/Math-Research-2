import sympy as s,json
K,z=s.symbols('K z',positive=True)
b=s.Rational(5,2)
crit=s.simplify(z**2*K**(-b)*s.sqrt(s.pi)*K**(b/2+s.Rational(1,4))*(z*z/2)**(s.Rational(1,4)-b/2))
assert s.simplify(crit-2*s.sqrt(s.pi)/K)==0
print(json.dumps(dict(sympy=s.__version__,critical_flank_factor=str(crit))))
