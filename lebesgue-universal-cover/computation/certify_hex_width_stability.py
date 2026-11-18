"""Bound the displacement to an exact hexagon-fitting copy in support distance."""
from fractions import Fraction as Q
import json
from pathlib import Path
import platform
import flint
from flint import arb, ctx

ctx.prec = 128
pi = arb.pi()
models = {}
for n in (3,5,7,9):
    a = pi/(2*n)
    r = 1/(2*a.cos())
    x = a if n % 3 == 0 else a/3
    amp = 3*r if n % 3 == 0 else r*(2*(2*a/3).cos()-1)
    models[f"B{n}"] = (amp*(1-x.cos())/x, r*a.sin())
for name,k,a5,b in (("F49",3,-Q(19,1200),Q(1,80)),("F50",15,Q(1,200),Q(3,1792))):
    def enclose(z):
        return arb(z.numerator)/z.denominator
    models[name] = (6*k*enclose(abs(b))/pi, enclose(5*abs(a5)+k*abs(b)))
rows = {}
for name,(slope,lip) in models.items():
    assert slope > 0
    constant = arb(4)/3 + 7*lip/slope
    assert constant < 45
    rows[name] = {"residual_per_radian_lower_enclosure":str(slope),
                  "support_lipschitz_enclosure":str(lip), "constant_enclosure":str(constant)}
base = Q(json.loads(Path(__file__).with_name("053-hex-fourier.json").read_bytes())["lower"])
cases = []
for delta,limit in [(Q(9,10**6),Q(".83674")), (Q(1,10**6),Q(".8379481"))]:
    lower = base/(1+90*delta)**2
    assert lower > limit
    cases.append({"width_excess":str(delta),"area_lower":str(lower),"strict_over":str(limit)})
print(json.dumps({"python":platform.python_version(),"python_flint":flint.__version__,
                  "precision":ctx.prec,"constant_bound":45,"models":rows,"cases":cases},indent=2))
