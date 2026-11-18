"""Enclose seven-body placements in irregular exact-width circumscribed hexagons."""
import json
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
import platform

import flint
from flint import arb, ctx

from certify_nine import FOURIER_TERMS, append_point, ar, arc_body, arithmetic, fourier_body, hull

K=2048
RHO_PAD=F(1,10**13)
NS=(3,5,7,9)


def support(j,alpha,pi):
    if j >= 4:
        return ar(F(1,2))+sum(ar(a)*(n*alpha).cos() for n,a in FOURIER_TERMS[j-4])
    n=NS[j]
    period=2*pi/n
    multiple=round(float(alpha.mid())/float(period.mid()))
    diff=abs(alpha-multiple*period)
    assert diff<pi/n
    half_vertex=pi/(2*n)
    radius=1/(2*half_vertex.cos())
    vertex=radius*diff.cos()
    arc=1-radius*(pi/n-diff).cos()
    if diff<half_vertex:
        return vertex
    if diff>half_vertex:
        return arc
    return arb.union(vertex,arc)


def determinant(u,v):
    return u[0]*v[1]-u[1]*v[0]


def verify(spec,pi):
    d1,d2=map(ar,spec["angle_perturbations"])
    theta=(arb(0),2*pi/3+d1,4*pi/3+d2)
    normals=[(a.cos(),a.sin()) for a in theta]
    w=[determinant(normals[(j+1)%3],normals[(j+2)%3]) for j in range(3)]
    assert all(t>0 for t in w)
    assert normals[1][1]>0

    points=[]
    for j in range(8*K):
        a=2*pi*j/(8*K)
        append_point(points,a.cos()/2,a.sin()/2)
    boxes=[]
    for j,c in enumerate(map(F,spec["rotation_root_proposals"])):
        def equation(phi):
            return sum(z*(support(j,a-phi,pi)-ar(F(1,2))) for z,a in zip(w,theta,strict=True))
        left,right=c-RHO_PAD,c+RHO_PAD
        fl,fr=equation(ar(left)),equation(ar(right))
        assert (fl<0 and fr>0) or (fl>0 and fr<0)
        phi=arb.union(ar(left),ar(right))
        tx=ar(F(1,2))-support(j,-phi,pi)
        ty=(ar(F(1,2))-support(j,theta[1]-phi,pi)-tx*normals[1][0])/normals[1][1]
        # At the root the third relation is forced by det-weight compatibility.
        if j<4:
            arc_body(points,K,pi,NS[j],phi,tx,ty)
        else:
            size=(8,16)[j-4]*K
            fourier_body(points,size,pi,FOURIER_TERMS[j-4],phi,tx,ty)
        boxes.append({
            "angle_root_bracket":[str(left),str(right)],
            "endpoint_signs":["+" if fl>0 else "-","+" if fr>0 else "-"],
            "tx":str(tx),
            "ty":str(ty),
        })
    poly=hull(points)
    res=arithmetic(poly,K)
    assert F(res["upper"]) < F("0.83809903")
    return {
        "angle_perturbations":spec["angle_perturbations"],
        "all_three_strip_widths_exactly_one":True,
        "certified_placement_boxes":boxes,
        "sample_count":len(points),
        "grid_count":len(poly),
        "area":res,
    }


def main():
    ctx.prec=160
    raw=Path(__file__).with_name("irregular_inputs.json").read_bytes()
    print(json.dumps({
        "python":platform.python_version(), "python_flint":flint.__version__,
        "resolution":K, "bits":ctx.prec, "source_sha256":sha256(raw).hexdigest(),
        "outputs":[verify(s,arb.pi()) for s in json.loads(raw)],
    },indent=2))


if __name__=="__main__":
    main()
