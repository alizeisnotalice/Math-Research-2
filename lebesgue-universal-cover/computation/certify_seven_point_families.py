"""Enclose a fixed seven-input hull using true arc samples and a rational parallel-body bound."""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import isqrt
from pathlib import Path

from certify_nine import S, append_point, ar, arb, arc_body, ctx, hull

HERE = Path(__file__).resolve().parent


def verify(k):
    ctx.prec = 160
    names = ["seven-point-families-pose.json","finite-point-leaf-input.json","multi-host-points.json"]
    raw = [(HERE/name).read_bytes() for name in names]
    params, p30, q11 = map(json.loads,raw)
    numeric = list(map(F,params["pose"]))
    assert len(numeric)==18 and numeric[-1] in (-1,1)
    p=list(map(ar,numeric))
    pi=arb.pi()
    pts=[]
    for i in range(8*k):
        a=pi*i/(4*k)
        append_point(pts,a.cos()/2,a.sin()/2)
    arc_body(pts,k,pi,3,arb(0),p[0],p[1])
    for n,j in [(5,2),(7,5),(9,8)]:
        arc_body(pts,k,pi,n,*p[j:j+3])
    for cloud,den_key,j in [(p30,"denominator",11),(q11,"point_denominator",14)]:
        xy=[(F(x,cloud[den_key]),F(y,cloud[den_key])) for x,y in cloud["points"]]
        assert len(xy)==(30 if j==11 else 11)
        assert max((x-a)**2+(y-b)**2 for x,y in xy for a,b in xy)<=1
        if j==14:
            cx=sum(x for x,y in xy)/len(xy)
            cy=sum(y for x,y in xy)/len(xy)
            xy=[(x-cx,numeric[-1]*(y-cy)) for x,y in xy]
        sn,cs=p[j].sin_cos()
        for x,y in xy:
            append_point(pts,cs*ar(x)-sn*ar(y)+p[j+1],sn*ar(x)+cs*ar(y)+p[j+2])
    vertices=hull(pts)
    edges=zip(vertices,vertices[1:]+vertices[:1])
    twice,perimeter=0,0
    for (x,y),(a,b) in edges:
        twice+=x*b-y*a
        perimeter+=isqrt((x-a)**2+(y-b)**2)+1
    area=F(twice,2*S*S)
    pi_up=F(22,7)
    rounding=F(2,S)
    epsilon=rounding+pi_up**2/(72*k*k)
    upper=area+F(perimeter,S)*epsilon+pi_up*epsilon**2
    # Rounded material lies within rounding of the true hull, which contains the half disk.
    lower=area/(1+2*rounding)**2
    assert 0<lower<upper<F(".837906")
    return dict(arc_samples=k,vertices=len(vertices),area_lower=str(lower),area_upper=str(upper),
                decimal_bracket=[f"0.{(lower * 10**16).__floor__():016d}",
                                 f"0.{(upper * 10**16).__ceil__():016d}"],precision_bits=ctx.prec,
                input_sha256={n:sha256(r).hexdigest() for n,r in zip(names,raw)})


if __name__=="__main__":
    print(json.dumps(verify(8192),indent=2))
