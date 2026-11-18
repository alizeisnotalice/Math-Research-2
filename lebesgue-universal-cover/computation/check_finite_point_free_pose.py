"""Bound the hull area at a specified free placement of the same six inputs."""
import argparse
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import json
import sys

from flint import arb,ctx

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from certify_nine import S, append_point, arc_body, ar
from certify_finite_point_leaf import finite_clouds,digest


def polygon(points):
    a=sorted(set(points))
    def chain(seq):
        out=[]
        for x,y in seq:
            while len(out)>=2:
                (u,v),(w,z)=out[-2:]
                if (w-u)*(y-v)-(z-v)*(x-u)>0:
                    break
                out.pop()
            out.append((x,y))
        return out
    return chain(a)[:-1]+chain(reversed(a))[:-1]


def bounds(material,error):
    hull=polygon(material)
    den=2*S*S
    area=F(sum(x*v-y*u for (x,y),(u,v) in zip(hull,hull[1:]+hull[:1])),den)
    perimeter=F(sum(isqrt((x-u)**2+(y-v)**2)+1 for (x,y),(u,v) in zip(hull,hull[1:]+hull[:1])),S)
    assert area>0 and max(x*x+y*y for x,y in hull)<S*S
    pi=F(22,7);a=F(2,S)
    low=area-2*pi*a-pi*a*a
    high=area+perimeter*error+pi*error*error
    return [str(low),str(high)],[float(low),float(high)]


def cloud_at(points,den,angle,tx,ty):
    out=[]
    s,c=angle.sin_cos()
    for x,y in points:
        append_point(out,tx+(c*x-s*y)/den,ty+(s*x+c*y)/den)
    return out


def verify(pose_path):
    ctx.prec=160
    data=json.loads((HERE/"finite-point-leaf-input.json").read_bytes())
    spec=json.loads(pose_path.read_bytes())
    vals=[ar(F(v)) for v in spec["pose"]]
    assert len(vals)==14
    geom=[(arb(0),arb(0),arb(0)),(arb(0),vals[0],vals[1])]
    geom.extend(tuple(vals[k:k+3]) for k in (2,5,8,11))
    assert all(abs(x)+abs(y)<ar(F(".2")) for _,x,y in geom)
    assert all(x*x+y*y<F(".7")**2*data["denominator"]**2 for x,y in data["points"])
    assert 1/(2*(arb.pi()/6).cos())<ar(F(".6"))
    P=cloud_at(data["points"],data["denominator"],*geom[-1])
    curved=list(P)
    density=10000
    for j,n in enumerate((0,3,5,7,9)):
        theta,tx,ty=geom[j]
        if n:
            arc_body(curved,(density+n-1)//n,arb.pi(),n,theta,tx,ty)
        else:
            for i in range(density):
                s,c=(2*arb.pi()*i/density).sin_cos()
                append_point(curved,c/2,s/2)
    pi=F(22,7)
    exact_c,decimal_c=bounds(curved,pi*pi/(4*density*density)+F(2,S))
    groups,_,_=finite_clouds()
    finite=list(P)
    for group,g in zip(groups,geom[:5],strict=True):
        finite.extend(cloud_at(group["points"],group["denominator"],*g))
    exact_f,decimal_f=bounds(finite,F(2,S))
    return dict(pose=spec["pose"],source=spec.get("source"),input_sha256=digest(HERE/"finite-point-leaf-input.json"),
                script_sha256=digest(Path(__file__)),curved_interval=exact_c,curved_decimal=decimal_c,
                finite_interval=exact_f,finite_decimal=decimal_f,arc_density=density)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("pose",type=Path)
    print(json.dumps(verify(p.parse_args().pose),indent=2))
