"""Exact inscribed polygon at one center distance and a dyadic volume ledger."""
import contextlib
import io
import json
import runpy
from fractions import Fraction as F
from pathlib import Path
import sys

base = Path(__file__).resolve().parent
with contextlib.redirect_stdout(io.StringIO()):
    union = runpy.run_path(str(base / "verify_projection_union.py"))
data = json.loads((base / "outer_domain_certificate.json").read_text())
scale = int(data["disk_polygon_denominator"])
disk = [tuple(map(int,p)) for p in data["disk_polygon_points"]]
assert all(4*(x*x+y*y)<=scale**2 for x,y in disk)
radius = F(283,600)
assert radius == F(1,2)-F(1,80)-F(19,1200)
distance = F(12,125)
# A single integer lattice represents both polygons exactly.
den = 300*scale
points = [(300*x,300*y) for x,y in disk]
shift = distance*den
assert shift.denominator == 1
points.extend((int(shift)+283*x,283*y) for x,y in disk)
ordered = sorted(set(points))
def cross(a,b,c):
    return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
def chain(seq):
    out = []
    for q in seq:
        while len(out)>=2 and cross(out[-2],out[-1],q)<=0:
            out.pop()
        out.append(q)
    return out[:-1]
h = chain(ordered)+chain(reversed(ordered))
area = F(sum(x*v-y*u for (x,y),(u,v) in zip(h,h[1:]+h[:1])),2*den**2)
assert area>F(419,500)

# Clip only Fourier cells lying wholly beyond the proven radius.
root = union["root"]
offset, high = root[2]
step, rem = divmod(high-offset,128)
assert rem==0
def near(lo,hi):
    return max(lo,-hi,0)
from functools import cache
@cache
def exterior_area(box):
    (x0,x1),(y0,y1)=box
    assert all((z-offset)%step==0 for z in (x0,x1,y0,y1))
    return step**2*sum((near(x,x+step)**2+near(y,y+step)**2)*distance.denominator**2
                       >= (scale*distance.numerator)**2
                       for x in range(x0,x1,step) for y in range(y0,y1,step))

extra=0
for leaf in union["current"]["leaves"]:
    if leaf["status"]=="retain":
        box=tuple(tuple(map(int,p)) for p in leaf["box_scaled"])
        b_remaining=union["volume"](box[:2])-union["prior_overlap"](box[:2])
        extra+=b_remaining*exterior_area(box[2:])
total=union["volume"](root)
remaining=F(union["remainder"]-extra,total)
assert 0<remaining<F(45377,524288)
print(json.dumps(dict(python=sys.version, radius=str(radius), distance=str(distance),
    polygon_area=str(area), added_fraction=str(F(extra,total)),
    remaining_fraction=str(remaining), excluded_fraction=str(1-remaining)),indent=2))
