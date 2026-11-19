""" independently pay every angular row using exact turns and direct endpoint rotations."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from math import isqrt
from pathlib import Path

from flint import arb,ctx
import numpy as np

N=384
U=2**22
Q=2**28
E=F(".000000994")
SIG=1-E
S=F(".999999038")


def read(p):return json.loads(p.read_bytes())
def dig(p):return sha256(p.read_bytes()).hexdigest()


def ball(x):
    x=F(x)
    return arb(x.numerator)/x.denominator


def lower(ball_value):
    v=ball_value.lower()
    return v


class CircleMap:
    def __init__(self, knots):
        assert len(knots)==N
        self.knots=knots
        ds=[]
        for j in range(N):
            a0,a1=knots[j],knots[(j+1)%N]
            v=1+N*ball(a1-a0)/(2*arb.pi())
            assert ball(".99")<v<ball("1.01")
            ds.append(v)
        self.speed=ds

    def shift(self, t):
        n0=N*t
        a0=n0.numerator//n0.denominator
        w=n0-a0
        idx=a0%N
        return (1-w)*self.knots[idx]+w*self.knots[(idx+1)%N]

    def receipts(self):
        for k in (384,1536):
            for j in range(k+1):
                t=F(j,k)
                assert self.shift(t)==self.shift(t+1)
                if k==1536 and j<k:
                    lo,jj=divmod(j,4)
                    d=F(jj,4)
                    expect=(1-d)*self.knots[lo]+d*self.knots[(lo+1)%N]
                    assert expect==self.shift(t)
        assert self.shift(F(0))==self.shift(F(1))
        assert self.shift(F(-1,1536))==self.shift(F(1535,1536))
        return dict(min_slope=min(z.lower() for z in self.speed).str(56),
                    max_slope=max(z.upper() for z in self.speed).str(56),
                    knot0=str(self.knots[0]), last_knot=str(self.knots[-1]),
                    largest_radian_shift=str(max(abs(b) for b in self.knots)))


def check_points(source, fit, actual, expected):
    raw=read(source)
    calibration=read(fit)
    unit=F(calibration["epsilon"])/calibration["den"]
    assert len(raw["points"])==len(actual["points"])==expected
    den=actual["denominator"]
    hand=([-1,1] if expected==11 else [1])
    for original, delta, final in zip(raw["points"],calibration["directions"],actual["points"],strict=True):
        assert all(F(x,raw["denominator"])+unit*y==F(z,den)
                   for x,y,z in zip(original,delta,final,strict=True))
    pts=actual["points"]
    d2=max(F((x-u)**2+(y-v)**2,den**2) for (x,y),(u,v) in combinations(pts,2))
    assert SIG**2<d2<S**2
    if expected==66:
        assert {tuple(p) for p in pts}=={(x,-y)for x,y in pts}
    else:
        reflected={tuple(p) for p in pts}
        assert reflected!={(x,-y) for x,y in pts}
    angles=calibration["angles"]
    assert len(angles)==len(hand)*N
    maps={eta:CircleMap([unit*v for v in angles[(eta==-1)*N:((eta==-1)+1)*N]])
          for eta in hand}
    return raw,maps,dict(diameter_squared=str(d2),strict_diameter_spare=str(S*S-d2),
                         maps={str(h):m.receipts() for h,m in maps.items()})


def candidate_at(xy, normal, phi, eta):
    x,y=xy
    vx,vy=normal
    sn,cs=phi
    return (vx*cs+vy*sn)*x+(vy*cs-vx*sn)*y*eta


def proof_for_bank(bank, old, spec, shape, actual, wall, maps):
    X=shape["points"]
    Y=actual["points"]
    den=actual["denominator"]
    dirs=bank["directions"]
    geometry=np.array(dirs)
    roots=1638 if spec in ["nativeQ","expandedQ"] else 5442
    assert wall["root_count"]==roots
    exact_d=[]
    for vx,vy in dirs:
        k=(vx*vx+vy*vy)*(Q//(2*U))**2
        n=isqrt(k)
        exact_d.append(n+int(4*k>(2*n+1)**2)+2-arb(k).sqrt())
        assert exact_d[-1]>0 and exact_d[-1]<3
    if spec=="strongQ":
        segments=[(N,N,1),(N,N,-1),(2*N,4*N,1),(2*N,4*N,-1)]
        assert len(bank["cells"])==18*N
    elif spec=="strongP":
        segments=[(0,N,1)]
        assert len(bank["cells"])==18*N
    elif spec=="nativeP":
        segments=[(0,N,1)]
        assert len(bank["cells"])==N
    else:
        segments=[(0,N,1),(0,N,-1)]
        assert len(bank["cells"])==N

    checks=0
    pressure=None
    curvature=None
    for off,K,eta in segments:
        cm=maps[eta]
        reference=np.array(X,dtype=np.float64)
        reference[:,1]*=eta
        for j in range(K):
            l,r=F(j,K),F(j+1,K)
            a0,a1=cm.shift(l),cm.shift(r)
            alpha_l=2*arb.pi()*ball(l)+ball(a0)
            alpha_r=2*arb.pi()*ball(r)+ball(a1)
            phi=(alpha_l.sin_cos(),alpha_r.sin_cos())
            k0=(N*l).numerator//(N*l).denominator
            speed=cm.speed[k0]
            assert abs((alpha_r-alpha_l)/(2*arb.pi()/K)-speed)<arb("1e-60")
            assert alpha_r-alpha_l<arb.pi()
            angle=2*np.pi*float((l+r)/2)
            cs,sn=np.cos(angle),np.sin(angle)
            axes=geometry@np.array([[cs,-sn],[sn,cs]])
            ids=(axes@reference.T).argmax(axis=1)
            read_cell=off+j
            cell=bank["cells"][read_cell]
            assert len(cell)==len(wall["rows"][read_cell])
            calls={}

            def projection(z, idx):
                if (z,idx) not in calls:
                    calls[z,idx]=[candidate_at(Y[idx],dirs[z],p,eta) for p in phi]
                return calls[z,idx]

            for ii,raw in enumerate(cell):
                keys,wt,g=raw[:3],raw[3:6],raw[6]
                assert min(wt)>0 and sum(wt)==g
                assert all(sum(w*dirs[k][t] for k,w in zip(keys,wt,strict=True))==0 for t in (0,1))
                comp=[projection(v,int(ids[v])) for v in keys]
                ev=[sum(w*c[i] for w,c in zip(wt,comp,strict=True))/(g*den*U) for i in (0,1)]
                if spec in ["nativeP","nativeQ","expandedQ"]:
                    c=old[j][ii][3:]
                    if eta<0:
                        c=[(-y)%11 for y in c]
                    alt=[projection(v,int(k)) for v,k in zip(keys,c,strict=True)]
                    alt_ev=[sum(w*k[i] for w,k in zip(wt,alt,strict=True))/(g*den*U) for i in (0,1)]
                    if min(alt_ev)>min(ev):
                        ev=alt_ev
                assert min(ev)>0
                shift=ball(E*F(U-2,U)*F(wall["rows"][read_cell][ii],wall["scale"]))
                if spec in ["strongP","strongQ"]:
                    shift+=ball(SIG)*sum(w*exact_d[k] for w,k in zip(wt,keys,strict=True))/(g*Q)
                qlim=([max(raw[7],z)for z in raw[8:10]] if len(raw)==10 else [raw[7]]*2)
                for q,h in zip(qlim,ev,strict=True):
                    v=(h+shift-ball(F(q,Q))).lower()
                    assert v>0
                    pressure=v if pressure is None else min(pressure,v)
                if len(raw)==10:
                    for x in ev:
                        t=(2*arb.pi()**2*speed**2*x-18*ball(F(min(raw[8:10]),Q))).lower()
                        assert t>0
                        curvature=t if curvature is None else min(curvature,t)
                checks+=1
    return dict(rows=checks, min_endpoint_payment=pressure.str(55),
                min_bow_margin=None if curvature is None else curvature.str(55),
                actual_segments=segments)


def main():
    p=ArgumentParser(description=__doc__)
    for k in ("comp","source","rows","output"):
        p.add_argument(k,type=Path)
    a=p.parse_args()
    ctx.prec=256
    f=read(a.source/"points.json")
    base={"P66":a.comp/"scaled-eight-50/p66.json","physicalQ":a.rows/"old/q.json"}
    groups={}
    rec={}
    for key,num,fit in [("P66",66,"P962.json"),("physicalQ",11,"Q962.json")]:
        groups[key]=check_points(base[key],a.source/fit,f[key],num)
        rec[key]=groups[key][2]
        rec[key]["map_sha256"]=dig(a.source/fit)
    specs=[
      ("strongP962","angular-66","grid","P66","strongP",9934),
      ("nativeP962","pleg","P","P66","nativeP",2493),
      ("strongQ962","angular-66","grid","physicalQ","strongQ",128850),
      ("nativeQ962","q-standard","Q","physicalQ","nativeQ",4516),
      ("expandedQ962","q-expanded","Qexp","physicalQ","expandedQ",7588),
    ]
    sums=read(a.source/"wall-sums.json")
    proofs={}
    for named,bank,wall,key,spec,count in specs:
        bank_path=a.rows/"old"/(bank+".json")
        wall_path=a.rows/f"walls-{wall}.json"
        b=read(bank_path)
        assert dig(wall_path)==sums[wall]
        w=read(wall_path)
        src,mp,r=groups[key]
        proposal=src.get("cells") if key=="P66" else b.get("proposals")
        info=proof_for_bank(b,proposal,spec,src,f[key],w,mp)
        assert count==info["rows"]==read(a.source/(named+".json"))["rows"]
        proofs[named]=info
        print(named,info,flush=True)
    out=dict(point_cloud_sha256=dig(a.source/"points.json"),directed_precision=ctx.prec,
             maps_and_clouds=rec,recomputed_rotation_rows=proofs)
    a.output.write_text(json.dumps(out,sort_keys=True,indent=2)+"\n")

if __name__=="__main__":
    main()
