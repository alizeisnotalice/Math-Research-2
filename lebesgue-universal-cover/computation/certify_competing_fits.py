"""Certify two small feasible hulls for the nine-mode test."""
import argparse
import json
from pathlib import Path
from fractions import Fraction as F
from math import isqrt

from flint import arb, ctx

from certify_nine import S, PI_UPPER, ar, arc_body, fourier_body, hull
from certify_five_harmonic import disk_points, placement_area
from certify_large_support import curvature, containment

def four_plus(center, terms, pi):
    k,count=2048,8192
    points=[]
    disk_points(points,8*k,pi,F(1,2))
    arc_body(points,k,pi,3,arb(0),*map(ar,center[:2]))
    for n,j in ((5,2),(7,5)):
        arc_body(points,k,pi,n,*map(ar,center[j:j+3]))
    fourier_body(points,count,pi,terms,*map(ar,center[8:]))
    p=hull(points)
    rows=list(zip(p,p[1:]+p[:1],strict=True))
    area=F(sum(a*d-b*c for (a,b),(c,d) in rows),2*S**2)
    bounds=[]
    for (a,b),(c,d) in rows:
        q=(a-c)**2+(b-d)**2
        z=isqrt(q)
        bounds.append(z+(z*z<q))
    perimeter=F(sum(bounds),S)
    eps=PI_UPPER**2*max(F(1,72*k*k),F(1,2*count*count))+F(2,S)
    return area+perimeter*eps+PI_UPPER*eps**2

def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument("experiment",type=Path)
    args=a.parse_args()
    ctx.prec=160
    pi=arb.pi()
    assert pi<ar(PI_UPPER)
    d=json.loads(args.experiment.read_bytes())
    terms=[(n,F(x)) for n,x in d["terms"]]
    assert curvature(terms,pi)==F(1,100)
    parameters=list(map(F,d["four_plus"]["params"]))
    assert len(parameters)==11
    upper4=four_plus(parameters,terms,pi)
    assert upper4<F(".836874")
    host=d["zero_gain_five"]
    width=containment({**host,"terms":d["terms"]},terms,pi)
    old=placement_area(list(map(F,host["params"])),pi,k=2048)
    assert F(old["upper"])<F(".837658")
    print(json.dumps({"four_plus_s_area_upper":str(upper4),"four_plus_s_area_display":float(upper4),
                      "five_with_s_strictly_inside":{"minimum_support_slack":width,
                            "area_upper":old["upper"],"area_display":old["display_bracket"]}},indent=2))

if __name__=="__main__":
    main()
