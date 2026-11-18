"""Enclose the rescaled five cores and normalize the diameter of the last group."""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import isqrt
from pathlib import Path
import sys

from certify_polygon_083731 import LAMBDA, S, ar, append_point, arb, ctx, lower_radius, diameter_squared

HERE = Path(__file__).resolve().parent
ctx.prec = 160
raw = (HERE / "polygon-083731-points.json").read_bytes()
report = json.loads((HERE / "polygon-083731-transfer.json").read_bytes())
assert sha256(raw).hexdigest() == report["points_sha256"] == "37446025cb4191da3ccd9a0283aa593eb6dea37560bcec4522de42f5f9a6320c"
groups = json.loads(raw)
diameter2 = diameter_squared(groups[-1])
unit = 10**18
t = F(isqrt(unit**2 * diameter2.denominator // diameter2.numerator), unit)
assert t*t*diameter2 <= 1 < (t+F(1, unit))**2*diameter2

pi = arb.pi()
for group in groups[:-1]:
    n, k = group["n"], group["samples_per_arc"]
    if n == 0:
        k += 1
        points = []
        for j in range(k):
            sine, cosine = (2*pi*j/k).sin_cos()
            append_point(points, cosine/2, sine/2)
        group.update(samples_per_arc=k, points=sorted(set(points)))
    half = pi/(2*n*k) if n else pi/k
    loss = (1-half.cos())*(1 if n else ar(F(1,2))) + ar(F(12,S))
    assert loss < ar(1-t*LAMBDA)*lower_radius(n,pi)
    assert diameter_squared(group) <= 1

last = groups[-1]
last["points"] = [[x*t.numerator, y*t.numerator] for x,y in last["points"]]
last["denominator"] *= t.denominator
assert diameter_squared(last) == t*t*diameter2
size = sum(len(group["points"]) for group in groups)
assert size == 8922
lower = t*t*F(report["area_floor"])
assert lower > F(".837310001265")
path = Path(sys.argv[1])
path.write_text(json.dumps(groups,separators=(",",":"))+"\n")
print(json.dumps(dict(points=size,scale=str(t),area_floor=str(lower),
                      points_sha256=sha256(path.read_bytes()).hexdigest()),indent=2))
