"""Check two triangular jet systems symbolically."""
import sympy as sp
x, t, b = sp.symbols("x t b", real=True)
f = t**2 * sp.exp(b*t**2/2) * sp.cosh(x*t)/2
rows = [sp.expand(sp.diff(f,t,2*j).subs(t,0)) for j in range(1,5)]
mat = sp.Matrix([[r.coeff(x,2*i) for i in range(4)] for r in rows])
assert mat.is_lower and mat.det() != 0
print("even zero-node Jacobian determinant =", mat.det())
C,S=sp.symbols("C S")
for j in range(5):
    differentiated = sp.diff(f,t,j).subs({sp.cosh(x*t): C, sp.sinh(x*t): S})
    top = sp.expand(differentiated).coeff(x,j).subs({C:sp.cosh(x*t), S:sp.sinh(x*t)})
    target= (sp.cosh(t*x) if j%2==0 else sp.sinh(t*x))*t**2*sp.exp(b*t**2/2)/2
    assert sp.simplify(top-target)==0
print("positive-node leading terms checked through order 4")
