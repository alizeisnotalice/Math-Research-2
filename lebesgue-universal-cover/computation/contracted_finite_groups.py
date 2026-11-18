"""Construct rational outer polygons after joint core and P30 contraction."""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import isqrt
from pathlib import Path
import sys

from certify_polygon_083731 import HERE, ORDERS, arb, ctx, diameter_squared, sampled_group

ctx.prec = 160
audit = json.loads((HERE / "p30-contraction-083731-audit.json").read_bytes())
source = (HERE / "finite-point-leaf-input.json").read_bytes()
assert sha256(source).hexdigest() == audit["P30_input_sha256"]
data = json.loads(source)
last = {"points": data["points"], "denominator": data["denominator"]}
d2 = diameter_squared(last)
assert d2 == F(audit["P30_diameter_squared_exact"])
precision = 10**18
scale = F(isqrt(precision**2 * d2.denominator // d2.numerator), precision)
factor = scale / (1-F(audit["P30_contraction"]))
core_scale = factor * F(audit["core_scale"])
groups = [sampled_group(n, arb.pi(), core_scale) for n in ORDERS]
last["denominator"] *= scale.denominator
last["points"] = [[x*scale.numerator, y*scale.numerator] for x,y in last["points"]]
groups.append(last)
diameters = [diameter_squared(group) for group in groups]
assert all(d <= 1 for d in diameters)
area = factor**2 * F(audit["integer_area_floor"])
assert area > F(".83731119527")
raw = (json.dumps(groups,separators=(",",":")) + "\n").encode()
Path(sys.argv[1]).write_bytes(raw)
print(json.dumps(dict(area=str(area), core_scale=str(core_scale),
                      group_sizes=[len(g["points"]) for g in groups],
                      squared_diameters=list(map(str,diameters)),
                      sha256=sha256(raw).hexdigest()),indent=2))
