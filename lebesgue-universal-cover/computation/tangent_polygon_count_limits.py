"""Integer lower bounds on vertices needed for a uniform dilation of diameter-one groups."""
import json
from fractions import Fraction as F
from flint import arb, ctx

ctx.prec = 160
lower, target = F('.83731'), F('.83730')
ratio = lower/target
dilation = (arb(ratio.numerator)/ratio.denominator).sqrt()
angle = (1/dilation).acos()
pi = arb.pi()
curves = {}
for n in (3,5,7,9):
    pieces = pi/(2*n*angle)
    k = int(float(pieces.mid())) + 1
    assert k-1 < pieces < k
    curves[str(n)] = n*(k-1)
# Perimeter inequalities force diam(P) >= cos(pi/(2m))/cos(pi/m).
def disk_dilation(m):
    return (pi/(2*m)).cos()/(pi/m).cos()

m = 3
while disk_dilation(m) > dilation:
    m += 1
assert disk_dilation(m) < dilation < disk_dilation(m-1)
print(json.dumps(dict(curved_lower=str(lower),target=str(target),dilation_limit=str(dilation),
                      disk_vertex_lower=m,arc_polygon_vertex_lower=curves,
                      five_group_vertex_lower=m+sum(curves.values()),
                      total_including_unchanged_thirty=m+sum(curves.values())+30),indent=2))
