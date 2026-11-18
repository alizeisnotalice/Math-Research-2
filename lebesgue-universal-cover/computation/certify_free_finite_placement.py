"""Rigorous area upper bound for one free placement of the five curved tests and rational points."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import isqrt
from pathlib import Path
import json
import platform

import flint
from flint import arb, ctx

from certify_nine import S, PI_UPPER, append_point, ar, arc_body, hull
from certify_five_harmonic import disk_points


def replay(input_data, placement, k):
    den=input_data['denominator']
    P=input_data['points']
    assert type(den) is int and den>0 and len(P)>=2
    assert all(len(p)==2 and all(type(z) is int for z in p) for p in P)
    d=max(sum((a-b)**2 for a,b in zip(p,t,strict=True)) for p,t in combinations(P,2))
    assert d<=den**2
    exact=[F(v) for v in placement]
    assert len(exact)==14 and k>=64
    xx=list(map(ar,exact))
    pi=arb.pi()
    assert pi<ar(PI_UPPER)
    coords=[]
    disk_points(coords,8*k,pi,F(1,2))
    arc_body(coords,k,pi,3,arb(0),*xx[:2])
    for n,j in ((5,2),(7,5),(9,8)):
        arc_body(coords,k,pi,n,*xx[j:j+3])
    cs,sn=xx[11].cos(),xx[11].sin()
    for a,b in P:
        x,y=ar(F(a,den)),ar(F(b,den))
        append_point(coords,cs*x-sn*y+xx[12],sn*x+cs*y+xx[13])
    polygon=hull(coords)
    ed=list(zip(polygon,polygon[1:]+polygon[:1],strict=True))
    area=F(sum(x*z-y*w for (x,y),(w,z) in ed),2*S*S)
    def ceil_distance(a,b):
        zz=sum((x-y)**2 for x,y in zip(a,b))
        m=isqrt(zz)
        return m+(m*m<zz)
    perimeter=F(sum(ceil_distance(a,b) for a,b in ed),S)
    eps=PI_UPPER**2/(72*k*k)+F(2,S)
    upper=area+perimeter*eps+PI_UPPER*eps**2
    return dict(point_count=len(P),point_max_pair_square_num=d,point_denom=den,
                exact_pose=list(map(str,exact)),samples=len(coords),arc_subdivisions=k,
                polygon_hull_vertices=len(polygon),rational_polygon_area=str(area),
                perimeter_upper=str(perimeter),parallel_padding=str(eps),strict_upper=str(upper),
                strict_upper_display=float(upper))


def main():
    ap=ArgumentParser(description=__doc__)
    ap.add_argument('input',type=Path)
    ap.add_argument('pose',type=Path,help='JSON array of 14 decimals or the numerical search result')
    ap.add_argument('--degree',type=int,default=4096)
    ap.add_argument('--target',type=F,default=F('.83730'))
    args=ap.parse_args()
    ctx.prec=160
    raw=args.input.read_bytes()
    poses=args.pose.read_bytes()
    data=json.loads(poses,parse_float=str)
    record=data['records'][0]['pose'] if isinstance(data,dict) else data
    result=replay(json.loads(raw),record,args.degree)
    result.update(target=str(args.target),below_target=F(result['strict_upper'])<args.target,
                  input_sha256=sha256(raw).hexdigest(),pose_source_sha256=sha256(poses).hexdigest(),
                  script_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
                  python=platform.python_version(),flint=flint.__version__,precision=ctx.prec)
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
