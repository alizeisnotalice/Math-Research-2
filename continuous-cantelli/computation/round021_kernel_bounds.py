"""Exact Gaussian envelope algebra for the exterior and shifted-derivative estimates."""
import sympy as s

u, w, z = s.symbols('u w z', positive=True)
p = s.exp(-u**2/(2*w))/s.sqrt(2*s.pi*w)
assert s.simplify(s.diff(p,w) - s.diff(p,u,2)/2) == 0
for du,dw,power,expected in (
    (1,0,2,-z**3), (0,1,3,(z**2-1)*z**3/2)
):
    expression = s.diff(p,u,du,w,dw).subs(w,u**2/z**2)
    normalized = expression*u**power*s.sqrt(2*s.pi)*s.exp(z**2/2)
    assert s.simplify(normalized-expected)==0
    print(s.expand(expected))
print("SymPy", s.__version__)
