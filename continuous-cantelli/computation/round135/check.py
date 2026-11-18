"""Verify the normalized increasing-kernel derivative polynomials exactly."""
import json
import sympy as s
w,v=s.symbols('w v')
p=s.Integer(1)
rows=[]
for j in range(15):
    assert p.subs(v,0)==w**j
    rows.append(str(p))
    # w=u/z and v=z^-2, while the unmarked kernel is exp(zu)/z.
    p=s.expand((w-v)*p-w*v*s.diff(p,w)-2*v*v*s.diff(p,v))
print(json.dumps(dict(sympy=s.__version__,polynomials=rows),indent=2))
