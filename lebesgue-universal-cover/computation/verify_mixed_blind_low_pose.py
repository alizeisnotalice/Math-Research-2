"""Certify that a five-body sublevel placement actually exists inside the prespecified blind box."""
from argparse import ArgumentParser
from hashlib import sha256
from fractions import Fraction as F
from pathlib import Path
import json
import sys
import numpy as np
from flint import ctx
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from audit_mixed_support_blind_box import q,det,dot,vertices,support_upper,ceil_rat,sh
from tomar_leaf_box import leaf_box

def polar_hull(vertices):
    points=sorted(vertices)
    def orient(a,b,c):
        return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
    def half(ordered):
        result=[]
        for x in ordered:
            while len(result)>=2 and orient(result[-2],result[-1],x)<=0:
                result.pop()
            result.append(x)
        return result[:-1]
    return half(points)+half(reversed(points))

def main():
    p=ArgumentParser(description=__doc__);p.add_argument('cert',type=Path);p.add_argument('pose',type=Path);a=p.parse_args()
    if not __debug__:raise SystemExit("Assertions required.")
    ctx.prec=192;z=json.loads(a.pose.read_bytes())
    assert [z['part'],z['pos'],z['coordinate_denominator']]==[1202,32259,2**44]
    manifest=(a.cert/'SHA256SUMS').read_bytes()
    assert sha256(manifest).hexdigest()=='56cbe5dfdadc6b5aa07cf6ec9acf4f637a3c4fa78db54df0ae2d78916d292410'
    pins={n.removeprefix('./'):h for h,n in (s.split() for s in manifest.decode().splitlines())}
    for f in ('part1202.json','part1202.npy'):assert sh(a.cert/f)==pins[f]
    code=np.load(a.cert/'part1202.npy',allow_pickle=False);prefix=np.r_[1,1+np.cumsum(np.where(code<11,1,-1))]
    c,w=leaf_box(prefix,code,z['part'],z['pos'],json.loads((a.cert/'part1202.json').read_bytes()))
    pose=[F(i,2**44) for i in z['pose']]
    assert all(F(x)-F(y)<=v<=F(x)+F(y) for x,y,v in zip(c,w,pose))
    V=[vertices(n,q(pose[3*i-1]) if i else q(0)) for i,n in enumerate((3,5,7,9))]
    offsets={};div=8192
    for k in range(div):
        x,y=div*div-k*k,2*div*k
        for _ in range(4):
            n=(x,y);ell=q(dot(n,n)).sqrt()
            coeff=[ell.upper()/2]
            for i,v in enumerate(V):
                coeff.append((q(dot(pose[3*i:3*i+2],n))+support_upper(v,n,ell)).upper())
            offsets[n]=ceil_rat(max(coeff))
            x,y=-y,x
    dual={(F(n[0],b),F(n[1],b)):n for n,b in offsets.items()}
    h=polar_hull(list(dual));normals=[dual[x]for x in h]
    outer=[]
    for i,n in enumerate(normals):
        m=normals[(i+1)%len(normals)];b=offsets[n];d=offsets[m];v=det(n,m);assert v>0
        outer.append(((b*m[1]-n[1]*d)/v,(n[0]*d-b*m[0])/v))
    for i,n in enumerate(normals):
        v=outer[i];u=outer[i-1]
        assert det(n,(v[0]-u[0],v[1]-u[1]))>0 and dot(v,n)==dot(u,n)==offsets[n]
    area=sum(det(v,outer[(i+1)%len(outer)]) for i,v in enumerate(outer))/2
    ub=F((area*10**12).__ceil__(),10**12)
    assert 0<ub<F('.83741'),float(ub)
    print(json.dumps(dict(part=1202,pos=32259,orientations=4*div,active_facets=len(normals),
        five_body_area_upper=str(ub),decimal_upper=f'{float(ub):.12f}',margin_below_083741=str(F('.83741')-ub),
        pose_sha256=sh(a.pose),source_sha256=sh(Path(__file__))),indent=2))
if __name__=='__main__':main()
