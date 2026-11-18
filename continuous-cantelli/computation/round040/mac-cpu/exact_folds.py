"""Exact folded-saddle algebra for a smooth divergent Stieltjes reference."""
import json
import platform
from pathlib import Path
import sympy as s

z,r=s.symbols("z r")
# Two denominator inequalities used uniformly over positive Laplace variables.
low_one=s.expand(1-(1-2*z)*(1+z)**2)
low_curvature=s.expand(1-3*z-(1-6*z)*(1+z)**3)
assert all(c>=0 for c in s.Poly(low_one,z).all_coeffs())
assert all(c>=0 for c in s.Poly(low_curvature,z).all_coeffs())
assert s.Rational(6,5)*(1-2*s.Rational(9,200))>1
assert 1-6*s.Rational(9,200)>0
x=s.Rational(1,2)+r/4+3*r*r/16
p=x-x*x*(1-r*x*x+2*r*r*x**4)
assert s.series(2*x*(1-2*r*x*x+6*r*r*x**4),r,0,3)==1+s.O(r**3)
assert s.series(-s.log(4*p),r,0,3)==-r/4-3*r**2/32+s.O(r**3)
out=dict(uniform_r_range=["0","1/8"],scaled_fold_bracket=["1/2","3/5"],
         uniform_endpoint_action=str(s.Rational(6,5)*(1-2*s.Rational(9,200))),
         uniform_curvature_floor=str(2*(1-6*s.Rational(9,200))),
         fold_expansion="x=1/2+r/4+3r^2/16+O(r^3)",
         normalized_fold_inverse_log="−r/4−3r^2/32+O(r^3)",
         python=platform.python_version(),sympy=s.__version__)
Path('exact_folds.json').write_text(json.dumps(out,indent=2)+"\n")
print(json.dumps(out,indent=2))
