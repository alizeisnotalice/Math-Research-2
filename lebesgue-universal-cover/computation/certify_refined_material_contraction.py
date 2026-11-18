"""Tighter common material contraction for the R187 rational point certificate."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt
from pathlib import Path
import json
import numpy as np

from certify_polygon_083731 import ORDERS, Q, ar, arb, ctx, lower_radius, verified_table
from certify_tomar_six_high import SourceTree
from certify_shallow_cap_transfer import diameter_squared

SCALE = F(".999999006")
SHRINK = F(".000000757")


def read(p):
    return json.loads(p.read_bytes())


def check_material(tree,counts,points):
    # In table units the exact center differs from (whole + frac/M) by <1/M.
    # Expand the squared norm about whole to keep vector operations in int64.
    M = 32768
    target = (Q*SCALE)**2
    pi = arb.pi()
    offsets = np.r_[0,np.cumsum(counts)]
    peaks = []
    for i,n in enumerate(ORDERS):
        body = points[:,offsets[i]:offsets[i+1]].astype(np.int64)
        assert np.abs(body).max() < Q
        if n==0:
            assert np.all(body==body[0])
            peak = int(np.square(body[0]).sum(axis=1).max())
            assert 4*peak < target
        else:
            fixed = n==3
            if fixed:
                assert np.all(body==body[0])
            end = F(0) if fixed else F(tree.root["hi"][3*i-4])
            amplitude = Q*ar(SCALE)/(2*(pi/(2*n)).cos())
            peak = 0
            for j in range(1 if fixed else 1025):
                for k in range(n):
                    t = ar(end*j/1024)+2*pi*k/n
                    sy,cx = t.sin_cos()
                    c = [amplitude*cx,amplitude*sy]
                    whole = np.array([round(float(v.mid())) for v in c],dtype=np.int64)
                    frac = np.array([round(float(((v-int(w))*M).mid())) for v,w in zip(c,whole)],dtype=np.int64)
                    assert max(abs(whole)) < Q
                    assert all(abs(v-int(w))<1 and abs((v-int(w))*M-int(f))<1
                               for v,w,f in zip(c,whole,frac))
                    d = body[j]-whole
                    linear = 2*(np.abs(d)-d*frac).sum(axis=1)
                    upper = np.square(d).sum(axis=1) - ((-linear)//M) + 2
                    peak = max(peak,int(upper.max()))
            assert peak < target
        peaks.append(peak)
    return peaks


def verify(comp,tree,table32,table64):
    ctx.prec=160
    receipt=read(comp/"signed-shallow-cap-083742-replay.json")
    assert receipt["complete"] and len(receipt["high_rows"])==66
    assert all(F(row[2])>=F(".837425") for row in receipt["high_rows"])
    goal=F(receipt["integer_goal"])
    assert goal==F(120684940798483225,2**57)
    tree=SourceTree(tree)
    src=read(comp/"tomar-six-high-input.json")
    seed=read(comp/"tomar-083742-seed.json")
    assert seed["tree_manifest_sha256"]==sha256((tree.path/"SHA256SUMS").read_bytes()).hexdigest()
    peaks={}
    for n,path in ((32,table32),(64,table64)):
        counts,points,stamp=verified_table(path,tree,n,src["tables"][str(n)])
        peaks[str(n)]=check_material(tree,counts,points)
        if n==32:
            material_point=tuple(map(int,points[0,counts[0]+12]))
            assert material_point==(81827130,81827130)
            radial=read(comp/"five-body-radial-affine.json")
            assert stamp==radial["source_table_sha256"]
            packet=comp/"five-body-radial-material.npz"
            assert sha256(packet.read_bytes()).hexdigest()==radial["table_sha"]
            with np.load(packet,allow_pickle=False) as f:
                assert np.array_equal(f["points"],points[0]) and np.array_equal(f["counts"],counts)
    s=SCALE; rho=SHRINK; r=F(1,2); R=F(".703"); h=F(".01335"); pi=arb.pi()
    assert all(F(k["integer_radius"],2**44)<F(".069") and F(k["integer_area_lower"],2**57)>goal for k in seed["radii"])
    for j,n in enumerate(ORDERS[1:]):
        k=3*j
        dist=min(F(tree.root["hi"][k]),F(tree.root["hi"][k+1]),-F(tree.root["lo"][k]),-F(tree.root["lo"][k+1]))
        ins=ar(s)*lower_radius(n,pi)
        assert pi*(ar(s*r)**2+ins**2)/2+ar(dist)*(ar(s*r)+ins)>ar(goal)
    sr=ar(s*r); rr=ar(R)
    assert (pi-(sr/rr).acos())*sr**2+sr*(rr**2-sr**2).sqrt()>ar(goal)
    assert s*s*F(".837295")+s*r*h*h/R>goal and R*R<2*(s*r)**2
    cloud=read(comp/"finite-point-leaf-input.json")
    den=cloud["denominator"]
    radius=F(".528619")
    assert all((F(x,den)+F(".01007"))**2+F(y,den)**2<radius**2 for x,y in cloud["points"])
    m=rho*radius/(1-s); a=r-F(".069"); v=(1+s)*r*r+(1-s)*m*m
    slope=1+s-(2*m*r+v)/(4*r*r)
    gap=4*(a-m)*s*s-(1+s)*h-max(0,(v-2*m*(r+h))/(2*r+h))
    assert m<a<r and slope>0 and gap>0
    q=F(".6")*(1-s)
    loss=F(22,7)*(2*R*q+q*q)
    assert F(".837425")-loss>goal
    points=comp/"shallow-cap-points.json"
    assert sha256(points.read_bytes()).hexdigest()==receipt["coordinate_sha256"]
    sets=read(points)
    assert all(x["denominator"]==1000000229098 for x in sets[:5])
    core=F(10**12-5,1000000229098)
    diam=diameter_squared(cloud)
    factor6=F(isqrt(10**36*diam.denominator//diam.numerator),10**18)
    assert sets[5]["points"]==[[factor6.numerator*x,factor6.numerator*y] for x,y in cloud["points"]]
    assert sets[5]["denominator"]==factor6.denominator*den
    f=min(core/s,factor6/(1-rho))
    lower=f*f*goal
    assert lower>F(".83742126")
    # Boundary of this specific transfer method, with its frozen radius estimate.
    smin=F(".99999900555")
    center_distance=ar(smin)/(2*(pi/6).cos())
    sy,cx=(4*pi/3).sin_cos()
    x,y=(ar(F(k,Q)) for k in material_point)
    assert (x-center_distance*cx)**2+(y-center_distance*sy)**2 > ar(smin)**2
    hmin=F(".01325")
    mmax=F(".4032")
    assert F(".837295")+r*hmin*hmin/R<goal
    assert (1+smin)-(2*a*r+2*r*r+(1-smin)*a*a)/(4*r*r)>0
    assert -4*smin*smin+2<0
    envelope=4*(a-mmax)-(1+smin)*hmin-max(0,((1+smin)*r*r-2*mmax*(r+hmin))/(2*r+hmin))
    assert envelope<0
    max_rho=mmax*(1-smin)/radius
    method_limit=goal*(factor6/(1-max_rho))**2
    assert lower<method_limit<F(".83742127165")
    return dict(target=str(goal),curve_scale=str(s),point_shrink=str(rho),height=str(h),
                material_scale_necessary_lower=str(smin),method_upper=str(method_limit),
                point_factor=str(f),material_squared_maxima=peaks,shallow_margin=str(gap),
                exceptional_loss=str(loss),lower=str(lower),display=float(lower),
                coordinate_sha256=receipt["coordinate_sha256"])


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("computation","tree","table32","table64"):
        p.add_argument(name,type=Path)
    a=p.parse_args()
    print(json.dumps(verify(a.computation,a.tree,a.table32,a.table64),indent=2))
