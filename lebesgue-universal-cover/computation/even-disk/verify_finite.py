"""Compute the new 2010 disk inradius and seven-group lower end by exact edge arithmetic."""
from argparse import ArgumentParser
from fractions import Fraction as F
from gzip import decompress
from hashlib import sha256
import json
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/"directional-seven"/"tangent-points.zip"))
import certify as t
SC,S,G=F(".999999083"),F(".99999909"),F(120696470013529293,2**57)

p=ArgumentParser(description=__doc__)
for k in ("out",):
    p.add_argument(k,type=Path)
a=p.parse_args()
from zipfile import ZipFile
tiers=HERE.parent/"two-center-seven"/"finite-thresholds.json"
with ZipFile(HERE/"finite-2010.zip") as z:
    raw=decompress(z.read("expected/disk2010.json.gz"))
cloud=json.loads(raw)
points=cloud["points"]
den=cloud["denominator"]
assert len(set(map(tuple,points)))==2010
hull=t.convex_hull(points)
d2=t.axis_max(points)[0]
assert len(hull)==2010 and d2==t.rotating_diameter2(hull)<den**2
factors=[]
for p,q in zip(hull,hull[1:]+hull[:1],strict=True):
    x,y=t.clockwise_edge_normal(p,q)
    dist=x*p[0]+y*p[1]
    assert dist>0
    factors.append(F(4*dist**2,den**2*(x*x+y*y)))
i=min(range(len(factors)),key=factors.__getitem__)
best=factors[i]
base=json.loads(tiers.read_bytes())["tiers"][".83750"]
curves=[F(r["included_scale"])**2 for r in base["curve_certificates"][1:]]
group_lower=min(G/S**2,G/SC**2*min([best]+curves))
assert group_lower>F(67,80)
seven=json.loads(decompress((HERE/"seven-2010.json.gz").read_bytes()))
older=json.loads(decompress((HERE.parent/"two-center-seven"/"tangent-910-points.json.gz").read_bytes()))
assert len(seven)==len(older)==7 and seven[0]["points"]==points and seven[0]["denominator"]==den and seven[1:]==older[1:]
out=dict(disk_source_coordinate_sha256=sha256(raw).hexdigest(),
         original_non_disk_certificate_sha256=sha256(tiers.read_bytes()).hexdigest(),
         disk_exact_max_squared_diameter=str(F(d2,den**2)),disk_exact_squared_inclusion=str(best),
         disk_limiting_hull_edge=i,new_counts=[2010]+base["actual_counts"][1:],
         joint_area_lower=str(group_lower),joint_area_lower_decimal=t.decimal_floor(group_lower))
a.out.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print(out["new_counts"],out["joint_area_lower_decimal"])
