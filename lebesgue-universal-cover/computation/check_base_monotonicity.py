"""Symbolic identity check; sign proof is in proofs/009-gibbs-base.md."""
import sympy as s
import json,sys
c=s.symbols('c',positive=True);f=s.sqrt(3)*(c/s.sqrt(3)-s.Rational(1,2))**2/(4*c*c-3)
d=(2*s.sqrt(3)*c-3)**2/(3*(4*c*c-3)**2)
assert s.simplify(s.diff(f,c)-d)==0
assert s.simplify(s.sqrt(3)/2-2*f.subs(c,1)-(2-2/s.sqrt(3)))==0
print(json.dumps({'python':sys.version,'sympy':s.__version__,'derivative':str(d),'zero_angle_identity':True,'status':'passed'},indent=2))
