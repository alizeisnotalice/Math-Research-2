"""Exclude Q11 on the R183 six-body escape host and certify a nearby area gap."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

from flint import arb,ctx

import certify_high_six_host_neighborhood as neighborhood
from certify_mixed_support_budget import Host,PI_UP,ar,cloud,det,dot,rotate
from certify_six_structural_escape import corners

BASE=F(".83742")
TARGET=F(".83753")
RADIUS=F(".00035")
PAD=F(".00226")
NORMAL_DEN=10**12
ANCHOR_SHA="fa4107fb115c05d559e801474457d7d985b25692b586088602930b4499f3d48b"

def digest(data):
    return sha256(data).hexdigest()

def read(path):
    return json.loads(path.read_bytes())

def normalize(seed):
    sign,k=seed
    assert type(sign) is int and sign in (-1,1)
    assert type(k) is int and abs(k)<=NORMAL_DEN
    t=F(k,NORMAL_DEN)
    u=(sign*(1-t*t)/(1+t*t),sign*2*t/(1+t*t))
    assert dot(u,u)==1
    return u

def weights(normals):
    if len(normals)==2:
        a,b=normals
        assert all(x+y==0 for x,y in zip(a,b,strict=True))
        return [F(1,2)]*2
    assert len(normals)==3
    a,b,c=normals
    q=det(b,c),det(c,a),det(a,b)
    assert all(x>0 for x in q) or all(x<0 for x in q)
    w=[x/sum(q) for x in q]
    assert all(x>0 for x in w) and sum(w)==1
    assert [sum(x*a[i] for x,a in zip(w,normals,strict=True)) for i in (0,1)]==[0,0]
    return w

class Padding:
    def __init__(self,pose,points,pi):
        self.source=pose
        self.points=points
        self.pi=pi
        self.memo={}

    def normal(self,seed):
        seed=tuple(seed)
        if seed in self.memo:
            return self.memo[seed]
        num=normalize(seed)
        u=tuple(map(ar,num))
        H=ar(F(1,2))
        for translation,vertices in self.source[0]:
            t=dot(translation,u)
            for c in vertices:
                value=dot(c,u)
                H=H.max((value+t).upper())
                if not value>ar(F(-1,2)):
                    H=H.max((1+value+t).upper())
        for p in self.source[1]:
            H=H.max(dot(p,u).upper())
        self.memo[seed]=(num,u,H)
        return self.memo[seed]

    def evaluate(self,left,right,hand,rows):
        mid=ar(left+right)*self.pi
        sin,cos=mid.sin_cos()
        half=ar(right-left)*self.pi
        assert 0<half<self.pi/2
        ns=[self.normal(r[:2]) for r in rows]
        w=weights([r[0] for r in ns])
        A=B=C=arb(0)
        for (_,_,i),(_,u,h),weight in zip(rows,ns,w,strict=True):
            assert type(i) is int and 0<=i<len(self.points)
            x,y=self.points[i]
            p=rotate((ar(x),ar(hand*y)),cos,sin)
            A+=ar(weight)*dot(p,u)
            B+=ar(weight)*det(p,u)
            C+=ar(weight)*h
        assert A>0
        return A*half.cos()-abs(B)*half.sin()-C

    def make(self,hand):
        import numpy as np
        from scipy.optimize import linprog

        def seed(a,b):
            sign=1 if a>=0 else -1
            return sign,round(NORMAL_DEN*sign*b/(1+sign*a))
        theta=np.arange(720)*2*np.pi/720
        seeds=[seed(a,b) for a,b in zip(np.cos(theta),np.sin(theta),strict=True)]
        tab=[self.normal(s) for s in seeds]
        normals=np.asarray([v[0] for v in tab],dtype=float)
        h=np.asarray([float(v[2]) for v in tab])
        constraints=np.c_[normals,-np.ones(len(tab))]
        cloud=np.asarray(self.points,dtype=float)*[1,hand]
        choices=dict(A_ub=constraints,bounds=[(None,None)]*3,method="highs",
                     options={"primal_feasibility_tolerance":1e-9,"dual_feasibility_tolerance":1e-9})
        def cell(i,k):
            a,b=FQ(i,k),FQ(i+1,k)
            x=float(a+b)*np.pi
            c,s=np.cos(x),np.sin(x)
            P=cloud@np.array([[c,s],[-s,c]])
            M=normals@P.T
            ys=M.argmax(axis=1)
            rhs=h-M[np.arange(len(ys)),ys]
            for objective in ([0,0,1],[2e-8,1e-8,1]):
                opt=linprog(objective,b_ub=rhs,**choices)
                assert opt.success
                ids=np.flatnonzero(-opt.ineqlin.marginals>1e-9)
                if len(ids) not in (2,3):
                    continue
                chosen=[[*seeds[j],int(ys[j])] for j in ids]
                try:
                    lower=self.evaluate(a,b,hand,chosen)
                except AssertionError:
                    continue
                if lower>ar(PAD):
                    yield [i,k,chosen]
                    return
            assert k<18,(i,k,hand,float(opt.fun))
            yield from cell(2*i,k+1)
            yield from cell(2*i+1,k+1)
        return [r for i in range(64) for r in cell(i,6)]

    def check(self,hand,data):
        end=F(0)
        floor=arb(1)
        for i,k,rows in data:
            a,b=FQ(i,k),FQ(i+1,k)
            assert a==end<b<=1
            value=self.evaluate(a,b,hand,rows)
            assert value>ar(PAD)
            end=b
            floor=floor.min(value)
        assert end==1
        return {"closed_cells":len(data),"certified_lower":str(floor),"display":float(floor)}

def FQ(a,k):
    return F(a,2**k)

def verify_p30(host,placed):
    floor=arb(1)
    for a,b in zip(host.material,host.material[1:]+host.material[:1],strict=True):
        u,v=b[1]-a[1],a[0]-b[0]
        c=ar(F(u*a[0]+v*a[1],host.den))
        length=arb(u*u+v*v).sqrt()
        for x,y in placed:
            floor=floor.min((c-u*x-v*y)/length)
    assert floor>ar(F(".00040"))
    return floor

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("computation",type=Path)
    p.add_argument("witness",type=Path)
    p.add_argument("--generate",action="store_true")
    args=p.parse_args()
    ctx.prec=160
    pi=arb.pi()
    assert pi<ar(PI_UP)
    base=args.computation
    files={k:base/v for k,v in {"anchor":"structural-six-escape-input.json",
                               "P30":"finite-point-leaf-input.json","Q11":"multi-host-points.json"}.items()}
    raw={k:read(p) for k,p in files.items()}
    info={k:digest(p.read_bytes()) for k,p in files.items()}
    assert info["anchor"]==ANCHOR_SHA
    assert raw["anchor"]["sources"]["finite-point-leaf-input.json"]==info["P30"]
    position=raw["anchor"]["pose"]
    phase=raw["anchor"]["hand"]
    assert len(position)==14 and phase in (1,-1)
    source=corners(position,raw["P30"],phase,pi)
    host=Host(position[:11],pi)
    clearance=verify_p30(host,source[1])
    Q=raw["Q11"]
    assert len(Q["points"])==11
    pts,diameter=cloud(Q["points"],Q["point_denominator"])
    disk=ar(F(1,2)); R=ar(F(".703"))
    assert (pi-(disk/R).acos())*disk**2+disk*(R**2-disk**2).sqrt()>ar(TARGET)
    neighborhood.LOW=BASE
    capsules=neighborhood.Caps(host,pts,pi,TARGET,RADIUS)
    padding=Padding(source,pts,pi)
    config={"inputs":info,"normal_denominator":NORMAL_DEN,"inner_hull_sha256":host.hull_hash,
            "six_family_lower":str(BASE),"augmented_target":str(TARGET),"augmented_radius":str(RADIUS),
            "padding_target":str(PAD)}
    hands=(1,-1)
    if args.generate:
        data=dict(config=config,caps={str(s):capsules.make(s) for s in hands},
                              padding={str(s):padding.make(s) for s in hands})
        tmp=args.witness.with_suffix(".pending")
        tmp.write_text(json.dumps(data,separators=(",",":"))+"\n")
        tmp.replace(args.witness)
    else:
        data=read(args.witness)
    assert data["config"]==config
    out=dict(complete=True,config=config,inner_polygon_area=str(host.area),
             P30_inside_inner_polygon_certificate=str(clearance),Q11_max_pair_distance_squared_numerator=diameter,
             Q11_denominator=Q["point_denominator"],witness_sha256=digest(args.witness.read_bytes()),
             two_sided_augmented_area={str(s):capsules.check(s,data["caps"][str(s)]) for s in hands},
             full_turn_Q11_padding={str(s):padding.check(s,data["padding"][str(s)]) for s in hands})
    print(json.dumps(out,indent=2))

if __name__=="__main__":
    main()
