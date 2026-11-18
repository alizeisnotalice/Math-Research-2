import sympy as s
D,r,eta=s.symbols('D r eta',positive=True)
lam=(D+eta)/(D-r)
rate=(D+eta)**2/(2*lam)-(eta+r)**2/(2*(lam-1))
assert s.simplify(rate-(D-r)**2/2)==0
print(s.factor(rate))
