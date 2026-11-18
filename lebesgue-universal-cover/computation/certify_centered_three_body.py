"""Enclose the minimum for concentric D, B3 and B5 over all orientations."""
import json
from fractions import Fraction as F
from pathlib import Path
from math import isqrt
from certify_rotation_inballs import hull
from rational_trig import PI, trig, scale

base = Path(__file__).parent
source = json.loads((base/'outer_domain_certificate.json').read_text())
den = int(source['disk_polygon_denominator'])
disk = [tuple(map(int,p)) for p in source['disk_polygon_points']]
assert all(4*(x*x+y*y)<=den**2 for x,y in disk)
def edges(vertices):
    return zip(vertices,vertices[1:]+vertices[:1])
# Fixed disk polygon contains radius r and misses the disk by at most .00002.
r=F(24999,50000)
for (x,y),(u,v) in edges(hull(disk)):
    c=x*v-y*u
    assert c>0 and c*c>=r*r*den**2*((x-u)**2+(y-v)**2)
def vertices(n):
    lo,hi=trig(scale(PI,F(1,2*n)),'cos')
    points=[]
    for j in range(n):
        arg=scale(PI,F(2*j,n))
        coords=[]
        for kind in ('cos','sin'):
            a,b=trig(arg,kind); vals=[p/(2*q) for p in (a,b) for q in (lo,hi)]
            k=round(den*(min(vals)+max(vals))/2)
            assert all(abs(F(k,den)-p)<=F(1,den) for p in vals)
            coords.append(k)
        points.append(tuple(coords))
    return points
tri=vertices(3); pent=vertices(5)
# tan(theta/2) <= tan(pi/30) < .106; verified without division near zero.
arg=scale(PI,F(1,30))
assert trig(arg,'sin')[1]<F(106,1000)*trig(arg,'cos')[0]
cells=512
lower=[]; upper=[]
for i in range(cells):
    t=F(106*(2*i+1),2000*cells)
    co=(1-t*t)/(1+t*t);si=2*t/(1+t*t)
    rotated=[(round(co*x-si*y),round(si*x+co*y)) for x,y in pent]
    h=hull(disk+tri+rotated)
    area=F(sum(x*v-y*u for (x,y),(u,v) in edges(h)),2*den**2)
    # Rotation speed with respect to t <=2, pentagon radius <1.
    err=F(106,1000*cells)+F(4,den)
    lower.append(area*(1-err/r)**2)
    # At the center there is no angular penalty. Every center (even past
    # the fundamental interval) is a valid orientation for the upper bound.
    upper.append(area*(1+(F(1,2)-r+F(4,den))/r)**2)
print(json.dumps(dict(cells=cells,lower=str(min(lower)),upper=str(min(upper)),
    lower_display=float(min(lower)),upper_display=float(min(upper)),
    best_center_index=upper.index(min(upper))),indent=2))
