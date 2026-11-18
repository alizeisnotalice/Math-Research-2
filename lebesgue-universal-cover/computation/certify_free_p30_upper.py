"""Strict Arb and rational area bracket for one pose of D, B3, B5, B7, B9, P30."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import isqrt
from pathlib import Path
import platform

import flint
from flint import arb, ctx

from certify_nine import S, append_point, ar, arc_body, hull


def verify(pose, cloud, k, clearance):
    assert len(pose) == 14 and k >= 64
    x = [ar(v) for v in pose]
    ctx.prec = 160
    pi = arb.pi()
    pi_upper = F(22,7)
    assert pi < ar(pi_upper)
    samples = []
    for i in range(8*k):
        theta = 2*pi*i/(8*k)
        append_point(samples,theta.cos()/2,theta.sin()/2)
    arc_body(samples,k,pi,3,arb(0),*x[:2])
    for n,j in ((5,2),(7,5),(9,8)):
        arc_body(samples,k,pi,n,*x[j:j+3])
    old_samples = len(samples)
    den = cloud["denominator"]
    pts = cloud["points"]
    assert den > 0 and len(pts) == len(set(map(tuple,pts))) == 30
    diameter = max((p[0]-q[0])**2+(p[1]-q[1])**2 for p in pts for q in pts)
    assert diameter <= den**2
    c,s = x[11].cos(),x[11].sin()
    for aa,bb in pts:
        a,b = ar(F(aa,den)),ar(F(bb,den))
        append_point(samples,x[12]+c*a-s*b,x[13]+s*a+c*b)
    poly = hull(samples)
    edges = list(zip(poly,poly[1:]+poly[:1],strict=True))
    area = F(sum(a[0]*b[1]-a[1]*b[0] for a,b in edges),2*S*S)
    assert area>0
    perim_up = 0
    for a,b in edges:
        d2 = sum((ai-bi)**2 for ai,bi in zip(a,b,strict=True))
        r = isqrt(d2)
        perim_up += r+(r*r<d2)
    perim_up = F(perim_up,S)
    endpoint = F(2,S)
    # Chord deviation <= circle radius * (subarc angle)^2/8, for
    # the disk and the unit arcs of B3/B5/B7/B9 respectively.
    sagitta = max(F(1,256),F(1,72))*pi_upper**2/k**2
    eps = endpoint+sagitta
    radius = F(59,100)
    if not all(F(u*u+v*v,S*S)<(radius-eps)**2 for u,v in poly):
        radius = F(7,10)
        assert all(F(u*u+v*v,S*S)<(radius-eps)**2 for u,v in poly)
    old_hull = hull(samples[:old_samples])
    # Every rounded P30 point lies an explicitly bounded distance inside
    # every edge of the rounded five-body hull. Both rounding errors are paid.
    assert clearance >= 0
    required = (clearance+2*endpoint)*S
    for (a,b) in zip(old_hull,old_hull[1:]+old_hull[:1],strict=True):
        dx,dy = b[0]-a[0],b[1]-a[1]
        square = dx*dx+dy*dy
        r = isqrt(square)
        length_bound = r+(r*r<square)
        det_min = min(dx*(q[1]-a[1])-dy*(q[0]-a[0]) for q in samples[old_samples:])
        if clearance:
            assert det_min > length_bound*required
    lower = area - 2*pi_upper*radius*endpoint - pi_upper*endpoint**2
    upper = area + perim_up*eps + pi_upper*eps**2
    assert upper > lower
    return {
        "pose_order": ["t3x","t3y","rho5","t5x","t5y","rho7","t7x","t7y",
                       "rho9","t9x","t9y","rhoP","tPx","tPy"],
        "pose": pose, "point_count": len(pts), "point_diameter_squared": str(F(diameter,den**2)),
        "arc_subdivisions": k, "sample_count": len(samples), "hull_vertices": len(poly),
        "grid_hull_area": str(area), "grid_hull_perimeter_upper": str(perim_up), "epsilon": str(eps),
        "area_lower": str(lower), "area_upper": str(upper),
        "display_bracket": [float(lower),float(upper)],
        "P30_support_clearance_lower": str(clearance), "hull_radius_upper": str(radius),
        "certified_below": [str(t) for t in (F(".8374"),F(".83745"),F(".83750"),F(".837515"),F(".83752"),F(".83765")) if upper<t],
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("pose",type=Path)
    p.add_argument("cloud",type=Path)
    p.add_argument("--subdivisions",type=int,default=4096)
    p.add_argument("--clearance",type=F,default=F(0))
    a = p.parse_args()
    ctx.prec = 160
    raw = a.cloud.read_bytes()
    pose = a.pose.read_bytes()
    d = verify(json.loads(pose),json.loads(raw),a.subdivisions,a.clearance)
    d.update({"cloud_sha256":sha256(raw).hexdigest(),"pose_sha256":sha256(pose).hexdigest(),
              "python":platform.python_version(),"flint":flint.__version__,"precision":ctx.prec})
    print(json.dumps(d,indent=2))

if __name__=="__main__":
    main()
