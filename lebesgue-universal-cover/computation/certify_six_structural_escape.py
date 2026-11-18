"""Certify one six-body host outside every R166 reference neighbourhood."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import isqrt
from pathlib import Path

import numpy as np
from flint import arb,ctx

from certify_nine import PI_UPPER,S,append_point,arc_body,ar,hull

N = 8192
J = 256
GRID = 10**12

def load(path):
    return json.loads(path.read_bytes())

def square(a,b):
    return (a[0]-b[0])**2+(a[1]-b[1])**2

def corners(pose,p30,hand,pi):
    if p30 is not None:
        d=p30["denominator"]
        assert hand in (-1,1)
        assert max(square(a,b) for a in p30["points"] for b in p30["points"])<d*d
    p=list(map(ar,pose))
    bodies=[]
    for n,i in ((3,-1),(5,2),(7,5),(9,8)):
        rho=p[i] if i>=0 else arb(0)
        t=p[i+1:i+3]
        R=1/(2*(pi/(2*n)).cos())
        assert R+(t[0]**2+t[1]**2).sqrt()<ar(F(3,5))
        v=[]
        for j in range(n):
            s,c=(rho+2*pi*j/n).sin_cos()
            v.append((R*c,R*s))
        bodies.append((t,v))
    points=[]
    if p30 is not None:
        s,c=p[11].sin_cos()
        for i,j in p30["points"]:
            x,y=ar(F(i,d)),ar(F(j*hand,d))
            points.append((c*x-s*y+p[12],s*x+c*y+p[13]))
    return bodies,points

def width_grid(source,num,pi,lower):
    bodies,points=source
    def sup(u,v):
        ans=ar(F(1,2))
        half=ar(F(-1,2))
        def keep(val):
            nonlocal ans
            ans=ans.max(val.lower() if lower else val.upper())
        for (x,y),vertices in bodies:
            t=u*x+v*y
            for a,b in vertices:
                z=u*a+v*b
                keep(z+t)
                admissible=z<half if lower else not z>half
                if admissible:
                    keep(1+z+t)
        for x,y in points:
            keep(u*x+v*y)
        return ans.lower() if lower else ans.upper()
    out=[]
    for i in range(num):
        v,u=(pi*i/num).sin_cos()
        w=(sup(u,v)+sup(-u,-v))*GRID
        count=w.lower().floor() if lower else w.upper().ceil()
        out.append(int(count.unique_fmpz()))
    return np.asarray(out,dtype=np.int64)

def area(pose,points,pi):
    k=2048
    p=list(map(ar,pose))
    samples=[]
    for j in range(8*k):
        s,c=(2*pi*j/(8*k)).sin_cos()
        append_point(samples,c/2,s/2)
    for n,i in ((3,-1),(5,2),(7,5),(9,8)):
        rho=p[i] if i>=0 else arb(0)
        arc_body(samples,k,pi,n,rho,*p[i+1:i+3])
    for q in points:
        append_point(samples,*q)
    poly=hull(samples)
    edges=list(zip(poly,poly[1:]+poly[:1],strict=True))
    aq=F(sum(a[0]*b[1]-a[1]*b[0] for a,b in edges),2*S**2)
    pq=F(sum(isqrt(square(a,b)-1)+1 for a,b in edges),S)
    near=F(2,S)
    far=near+PI_UPPER**2/(72*k*k)
    rad=F(3,5)
    assert max(square((0,0),a) for a in poly)<(S*(rad-far))**2
    return (aq-2*PI_UPPER*rad*near-PI_UPPER*near**2,
            aq+pq*far+PI_UPPER*far**2),poly,far

def merge_indices(indices):
    groups=[]
    for j,a in enumerate(indices):
        a=int(a)
        if groups and groups[-1][-1]==a:
            groups[-1][1]=j+1
        else:
            groups.append([j,j+1,a])
    return groups

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("comp",type=Path)
    p.add_argument("input",type=Path)
    p.add_argument("--write-geometry",type=Path)
    a=p.parse_args()
    ctx.prec=160
    pi=arb.pi()
    assert pi<ar(PI_UPPER)
    spec=load(a.input)
    for name,sha in spec["sources"].items():
        assert sha256((a.comp/name).read_bytes()).hexdigest()==sha
    get=lambda n:load(a.comp/n)
    p30=get("finite-point-leaf-input.json")
    refs=[get("free-p30-reoptimized-pose.json")[:11],get("free-p30-placement.json")[:11],
          get("competing-fits-input.json")["zero_gain_five"]["params"]]
    x=spec["pose"]
    assert len(x)==14
    domain=corners(x,p30,spec["hand"],pi)
    fixed=width_grid(domain,J,pi,False)
    bracket,poly,error=area(x,domain[1],pi)
    assert F(".8375136706688")<bracket[0]<bracket[1]<F(".83751378")<F(".83753")

    angular=F(6,5)*PI_UPPER/(2*N)
    js=np.arange(J,dtype=np.int64)*(N//J)
    steps=np.arange(N,dtype=np.int64)
    tested=[]
    for label,ref,radius in zip(("new_P30","old_P30","escape108"),refs,map(F,spec["radii"]),strict=True):
        width=width_grid(corners(ref,None,1,pi),N,pi,True)
        hands={}
        for h in (1,-1):
            ids=h*(js[:,None]-steps[None,:])%N
            gap=width[ids]-fixed[:,None]
            ids=np.argmax(gap,axis=0)
            pixels=gap[ids,steps]
            lower=(F(int(pixels.min()),GRID)-angular)/2
            assert lower>radius
            hands[str(h)]={"radius_lower":str(lower),"display":float(lower),
                           "weakest_rotation_index":int(np.argmin(pixels)),
                           "witness_intervals":merge_indices(ids)}
        tested.append({"reference":label,"declared_radius":str(radius),"hands":hands})
    if a.write_geometry:
        proof={"pose":x,"hand":spec["hand"],"P30":dict(points=p30["points"],denominator=p30["denominator"]),
               "placed_P30_arb":[list(map(str,v)) for v in domain[1]],
               "reference_poses":refs,"polygon":poly,"denominator":S,
               "polygon_to_true_host_Hausdorff_upper":str(error)}
        a.write_geometry.write_text(json.dumps(proof,separators=(",",":"))+"\n")
    print(json.dumps({"complete":True,"input_sha256":sha256(a.input.read_bytes()).hexdigest(),
        "area_bracket":list(map(str,bracket)),"area_display":list(map(float,bracket)),
        "rotations_each_hand":N,"directions_per_host":J,"width_lipschitz_error":str(angular),
        "reference_noninclusions":tested},separators=(",",":")))

if __name__=="__main__":
    main()
