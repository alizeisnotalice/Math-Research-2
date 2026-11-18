"""Outward rational Green integral for the fixed combined AH cap."""
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
order, cells = 14, 24
original_args = sys.argv[:]
sys.argv = [str(coordinates), sigma, str(order), str(cells)]
geo = {"__file__": str(coordinates)}
exec(compile(prefix, str(coordinates), "exec"), geo)
sys.argv = original_args
I, Jet = geo["I"], geo["Jet"]
add, sub, mul, dot, cross, perp = (geo[k] for k in ("add", "sub", "mul", "dot", "cross", "perp"))
pair, trig = geo["pair"], geo["trig"]
z, D = geo["z"], geo["direction"]
pivot_parameter = F("0.38")
angle = F("0.00033")
W = add(z, mul(I(pivot_parameter), D))
length = dot(D, D).sqrt()
b = mul(-1 / length, D)
cosine, sine = (I(*trig.trig((angle, angle), key)) for key in ("cos", "sin"))
line_direction = add(mul(cosine, b), mul(sine, perp(b)))


def curve(position):
    return pair(add(z, mul(I(position), D)))


def curve_series(position):
    w = Jet(position)
    w.c[1] = I(1)
    return pair(add(tuple(map(Jet, z)), mul(w, tuple(map(Jet, D)))))


def cut(position):
    return cross(line_direction, sub(curve(position), W))


start, finish = F("0.037"), F("0.038")
root_check = [cut(start), cut(finish)]
assert root_check[0].hi < 0 < root_check[1].lo
iterations = 112
for _ in range(iterations):
    mid = (start + finish) / 2
    value = cut(mid)
    if value.hi < 0:
        start = mid
    else:
        assert value.lo > 0
        finish = mid

local_jets = curve_series(I(start, finish))
slope = cross(line_direction, tuple(p.c[1] for p in local_jets))
assert slope.lo > 0

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
width = start / cells
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

# The unknown remaining interval has a length between zero and finish-start.
end_density = density(I(start, finish))[0]
endpoint_correction = I(0, finish - start) * end_density
cap_area = (integral + endpoint_correction) / 2
assert cap_area.lo > F("0.00000004056255211811389")

baseline = here / "eh-cap-switch-360.json"
eh = json.loads(baseline.read_bytes())["summary"]
assert F(eh["sigma_degrees"]) == F(sigma) and F(eh["delta"]) == F("0.00036")
base = I(*map(F, eh["candidate_area_interval_exact"]))
body_area = base - cap_area


def decimal_enclosure(value, places):
    scale = 10**places
    lo = (scale * value.lo).__floor__()
    hi = (scale * value.hi).__ceil__()

    def decimal(n):
        return f"{n//scale}.{n % scale:0{places}d}"

    return [decimal(lo), decimal(hi)]


print(json.dumps({
    "python": sys.version,
    "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "coordinates_sha256": hashlib.sha256(source).hexdigest(),
    "EH_baseline_sha256": hashlib.sha256(baseline.read_bytes()).hexdigest(),
    "sigma_degrees_exact": sigma, "EH_depth_radians_exact": "0.00036",
    "AH_pivot_exact": str(pivot_parameter), "AH_angle_radians_exact": str(angle),
    "AH_crossing_parameter_exact_interval": [str(start), str(finish)],
    "AH_crossing_parameter_decimal_interval": decimal_enclosure(I(start, finish), 30),
    "AH_crossing_slope_lower_exact": str(slope.lo),
    "bisection_iterations": iterations, "Taylor_order": order, "Taylor_cells": cells,
    "Taylor_remainder_abs_upper_exact": str(remainder_abs),
    "endpoint_integral_correction_exact": endpoint_correction.out(),
    "AH_cap_area_exact_interval": cap_area.out(),
    "AH_cap_area_decimal_interval": decimal_enclosure(cap_area, 30),
    "combined_candidate_area_exact_interval": body_area.out(),
    "combined_candidate_area_decimal_interval": decimal_enclosure(body_area, 30),
    "passed": True,
}, indent=2))
