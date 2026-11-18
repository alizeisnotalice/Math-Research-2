"""Independent exact Hessian identity for the two Wick saddle signs."""
import json
from pathlib import Path
import platform
import sympy as s

R,U,C,E,nu,alpha,beta=s.symbols("R U C E nu alpha beta",positive=True)
J=s.Matrix([[C+E,-E],[-E,E]])
g=2*nu*alpha*s.Matrix([R,U])
factor=s.simplify((g.T*J.inv()*g)[0])
want=4*nu**2*alpha**2*((R+U)**2/C+U**2/E)
assert s.simplify(factor-want)==0
shift=-J.inv()*g
substituted=s.simplify(2*nu*beta*(R+U)**2+(g.T*shift)[0]+(shift.T*J*shift)[0]/2)
assert s.simplify(substituted-(2*nu*beta*(R+U)**2-want/2))==0
Path("saddle_symbolic.json").write_text(json.dumps(dict(
    hessian_determinant=str(J.det()),gradient_quadratic=str(s.factor(factor)),
    quadratic_phase=str(s.factor(substituted)),python=platform.python_version(),sympy=s.__version__),indent=2)+"\n")
print("exact Hessian verified")
