"""Exact maximum AH chord cap allowed by the compulsory boundary chain."""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys

here = Path(__file__).resolve().parent
coordinates = here / "certify_gibbs_area_refined.py"
source = coordinates.read_bytes()
prefix, _ = source.split(b"width = F(1,cells)")
sigma = "1.5497"
order, cells = 10, 64
original_args = sys.argv[:]
sys.argv = [str(coordinates), sigma, str(order), str(cells)]
geo = {"__file__": str(coordinates)}
exec(compile(prefix, str(coordinates), "exec"), geo)
sys.argv = original_args
I, Jet = geo["I"], geo["Jet"]
add, sub, mul, dot, cross, perp = (geo[k] for k in ("add", "sub", "mul", "dot", "cross", "perp"))
pair, trig = geo["pair"], geo["trig"]
z, D = geo["z"], geo["direction"]
critical=F("0.429")
W=add(z,mul(I(critical),D))

def curve_series(position):
    w=Jet(position)
    w.c[1]=I(1)
    return pair(add(tuple(map(Jet,z)),mul(w,tuple(map(Jet,D)))))

def density(position):
    x, y = curve_series(position)
    x.c[0] = x.c[0] - W[0]
    y.c[0] = y.c[0] - W[1]
    return [
        sum((j + 1) * (x.c[k-j] * y.c[j+1] - y.c[k-j] * x.c[j+1])
            for j in range(k + 1))
        for k in range(order + 1)
    ]


# Integrate from 0 to the rational lower endpoint. Every Taylor remainder
# coefficient is evaluated at the whole integration cell.
width = critical / cells
integral = I(0)
remainder_abs = F(0)
for cell in range(cells):
    a, b = cell * width, (cell + 1) * width
    coeff = density(I((a + b) / 2))
    remainder = density(I(a, b))[order]
    for k in range(0, order, 2):
        integral += (width ** (k+1) / ((k+1)*2**k)) * coeff[k]
    factor = width ** (order+1) / ((order+1)*2**order)
    integral += factor * remainder
    remainder_abs += factor * max(abs(remainder.lo), abs(remainder.hi))

cap_area=integral / 2
compulsory=here/"ah-compulsory-tail.json"
witness=json.loads(compulsory.read_bytes())
assert F(witness["sigma_degrees_exact"])==F(sigma)
assert F(witness["compulsory_tail_parameter_range_exact"][0])==critical
original=here/"gibbs-1p5497-refined-18-128.json"
eh_input=here/"eh-cap-switch-joint-552p1-revisit-ah.json"
old_area=I(*map(F,json.loads(original.read_bytes())["area_interval"]))
eh_area=I(*map(F,json.loads(eh_input.read_bytes())["summary"]["cut_gain_area_interval_exact"]))
area_floor=old_area-eh_area-cap_area
cap_lo,cap_hi=F("0.0000074358439921451580997320"),F("0.0000074358439921451580997321")
assert cap_lo<cap_area.lo<=cap_area.hi<cap_hi
assert area_floor.lo > F("0.8440860344")
print(json.dumps({
    "python":sys.version,
    "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "compulsory_chain_input_sha256":hashlib.sha256(compulsory.read_bytes()).hexdigest(),
    "EH_area_input_sha256":hashlib.sha256(eh_input.read_bytes()).hexdigest(),
    "Gibbs_area_input_sha256":hashlib.sha256(original.read_bytes()).hexdigest(),
    "compulsory_parameter_exact":str(critical),
    "Taylor_cells":cells,"Taylor_order":order,
    "maximum_AH_cap_area_interval_exact":cap_area.out(),
    "maximum_AH_cap_area_outward_decimal_interval":["0.0000074358439921451580997320",
                                                    "0.0000074358439921451580997321"],
    "restricted_body_area_floor_interval_exact":area_floor.out(),
    "restricted_body_area_floor_rounded_down":"0.8440860344",
    "passed":True,
},indent=2))
