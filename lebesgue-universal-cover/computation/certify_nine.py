"""Arb and integer certificates for fixed eight- and nine-body hull areas."""

import argparse
import json
import platform
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt
from pathlib import Path

import flint
from flint import arb, ctx

S = 10**12
PI_UPPER = F(22,7)
FOURIER_TERMS = [
    ((3,F(1,80)),(5,-F(19,1200))),
    ((5,F(1,200)),(15,F(3,1792))),
]
PI_QUALITY = 160


def ar(value):
    r = F(value)
    return arb(r.numerator) / r.denominator


def append_point(points,x,y):
    q = tuple(round(float(z.mid())*S) for z in (x,y))
    for value,coordinate in zip((x,y),q,strict=True):
        assert ar(F(coordinate-1,S)) < value < ar(F(coordinate+1,S))
    points.append(q)


def arc_body(points,k,pi,n,rho,tx,ty):
    radius = 1 / (2*(pi/(2*n)).cos())
    for ell in range(n):
        phi = rho+2*ell*pi/n
        cx,cy=(tx+radius*phi.cos(),ty+radius*phi.sin())
        for j in range(k+1):
            angle = phi + pi*(1 - ar(F(1,2*n)) + ar(F(j,n*k)))
            append_point(points,cx+angle.cos(),cy+angle.sin())


def polygon_body(points,pi,n,rho,tx,ty):
    radius = 1 / (2*(pi/(2*n)).cos())
    for j in range(n):
        theta = rho+2*j*pi/n
        append_point(points,tx+radius*theta.cos(),ty+radius*theta.sin())


def fourier_body(points,count,pi,coeffs,rho,tx,ty):
    for j in range(count):
        theta = 2*j*pi/count
        h=ar(F(1,2))+sum(ar(a)*(m*theta).cos() for m,a in coeffs)
        hp=-sum(ar(m*a)*(m*theta).sin() for m,a in coeffs)
        c,s=(theta+rho).cos(),(theta+rho).sin()
        append_point(points,tx+h*c-hp*s,ty+h*s+hp*c)


def cross(a,b,c):
    return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])


def hull(points):
    def half(seq):
        out = []
        for p in seq:
            while len(out)>1 and cross(out[-2],out[-1],p)<=0:
                out.pop()
            out.append(p)
        return out[:-1]
    ordered=sorted(set(points))
    return half(ordered)+half(reversed(ordered))


def arithmetic(poly,k):
    edges=list(zip(poly,poly[1:]+poly[:1],strict=True))
    twice=sum(a[0]*b[1]-a[1]*b[0] for a,b in edges)
    assert twice>0
    polygon_area=F(twice,2*S**2)
    perimeter_integers=0
    for p,q in edges:
        square=sum((a-b)**2 for a,b in zip(p,q,strict=True))
        length=isqrt(square)
        perimeter_integers+=length+(length*length<square)
    perimeter_upper=F(perimeter_integers,S)

    # The old Fourier curvature and derivative have respective sup bounds
    # .98 and 2.2, and the new values are .995 and 6.225.
    sup = (F(5,2),F(32,5))
    for coeffs,m in zip(FOURIER_TERMS,sup,strict=True):
        radial_wave = sum((n*n-1)*abs(a) for n,a in coeffs)
        radial_derivative = sum(n*(n*n-1)*abs(a) for n,a in coeffs)
        assert 0 < F(1,2)-radial_wave
        assert (F(1,2)+radial_wave)**2+radial_derivative**2 < m*m
    interpolation_constants=(
        F(1,2)/(8*4**2), F(1)/(8*3**2),
        sup[0]/(8*4**2), sup[1]/(8*8**2),
    )
    chord_error=max(interpolation_constants)*PI_UPPER**2/k**2
    assert chord_error == F(5,256)*PI_UPPER**2/k**2
    endpoint_error=F(2,S)
    epsilon=chord_error+endpoint_error

    # H subset C + endpoint_error B, and C subset H + epsilon B.
    # Bound the true hull's perimeter for the reverse Steiner comparison.
    radius=F(3,5)
    if not all(F(x*x+y*y,S*S) < (radius-epsilon)**2 for x,y in poly):
        radius=F(2)
        assert all(F(x*x+y*y,S*S) < (radius-epsilon)**2 for x,y in poly)
    assert radius > epsilon
    lower=polygon_area - 2*PI_UPPER*radius*endpoint_error - PI_UPPER*endpoint_error**2
    upper=polygon_area + perimeter_upper*epsilon + PI_UPPER*epsilon**2
    return {
        "grid_polygon_area":str(polygon_area),
        "true_hull_radius_upper":str(radius),
        "epsilon":str(epsilon),
        "hull_perimeter_upper":str(perimeter_upper),
        "lower":str(lower),
        "upper":str(upper),
        "display_bracket":[float(lower),float(upper)],
    }


def certify(params,k):
    exact = [ar(z) for z in params]
    pi = arb.pi()
    assert pi<ar(PI_UPPER)
    points=[]
    for j in range(8*k):
        angle=2*j*pi/(8*k)
        append_point(points,angle.cos()/2,angle.sin()/2)
    arc_body(points,k,pi,3,arb(0),*exact[:2])
    for n,j in ((5,2),(7,5),(9,8)):
        arc_body(points,k,pi,n,*exact[j:j+3])
    for n,j in ((11,11),(19,14)):
        polygon_body(points,pi,n,*exact[j:j+3])
    fourier_body(points,8*k,pi,FOURIER_TERMS[0],*exact[17:20])
    if len(exact)==23:
        fourier_body(points,16*k,pi,FOURIER_TERMS[1],*exact[20:23])
    p=hull(points)
    return {
        "params":params, "sample_count":len(points),"polygon_vertices":len(p),
        **arithmetic(p,k),
    },p


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source",type=Path)
    parser.add_argument("candidate",type=Path)
    parser.add_argument("output",type=Path)
    parser.add_argument("--degree",type=int,default=1024)
    args=parser.parse_args()
    ctx.prec=PI_QUALITY
    src_raw=args.source.read_bytes()
    opt_raw=args.candidate.read_bytes()
    source=json.loads(src_raw)
    candidate=json.loads(opt_raw)
    assert candidate["source_sha256"]==sha256(src_raw).hexdigest()
    base=source["parameter_decimals_interpreted_exactly"]
    joint=candidate["best"]["params"]
    fixed=base+candidate["best_fixed_eight"]["params"][20:]
    assert len(base)==20 and len(joint)==len(fixed)==23 and args.degree>=64
    output={}
    for label, params in (("old_eight",base),("old_eight_with_new_body",fixed),
                         ("reoptimized_nine",joint),("reoptimized_eight_subset",joint[:20])):
        out, poly=certify(params,args.degree)
        if label=="reoptimized_nine":
            out["hull_grid_points"]=poly
        output[label]=out
        print(label,out["display_bracket"],flush=True)
    def endpoint(label,end):
        return F(output[label][end])
    assert endpoint("reoptimized_nine","upper") < F(8383,10000)
    assert endpoint("reoptimized_nine","upper") < endpoint("old_eight","lower")
    assert endpoint("old_eight","upper") < endpoint("old_eight_with_new_body","lower")
    assert endpoint("reoptimized_eight_subset","upper") < endpoint("reoptimized_nine","lower")
    result={
        "degree":args.degree,
        "bits":PI_QUALITY,
        "python":platform.python_version(),
        "python_flint":flint.__version__,
        "grid_denominator":S,
        "script_sha256":sha256(Path(__file__).read_bytes()).hexdigest(),
        "old_certificate_sha256":sha256(src_raw).hexdigest(),
        "numerical_search_sha256":sha256(opt_raw).hexdigest(),
        "disjoint_comparisons":{
            "reoptimized_nine_cheaper_than_old_eight_by_at_least":str(
                endpoint("old_eight","lower")-endpoint("reoptimized_nine","upper")),
            "frozen_eight_costs_more_after_placing_ninth_by_at_least":str(
                endpoint("old_eight_with_new_body","lower")-endpoint("old_eight","upper")),
            "ninth_body_increases_area_at_reoptimized_placement_by_at_least":str(
                endpoint("reoptimized_nine","lower")-endpoint("reoptimized_eight_subset","upper")),
        },
        "certified_cases":output,
        "scope":"Exact-area brackets for specified placements and an upper bound for a finite nine-body test family.",
    }
    args.output.write_text(json.dumps(result,indent=2)+"\n")


if __name__=="__main__":
    main()
