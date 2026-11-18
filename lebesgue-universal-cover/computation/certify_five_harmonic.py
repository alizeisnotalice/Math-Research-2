"""Replay the five-body harmonic exclusion and its finite-family comparison."""
import argparse
import json
import platform
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt, lcm
from pathlib import Path

import flint
from flint import arb, ctx

import certify_nine
from certify_nine import S, PI_UPPER, append_point, ar, arc_body, arithmetic, cross, fourier_body, hull

# D, B3, B5, B7, B9, followed by the new smooth width-one body.
TERMS = ((3,F(1,100)),(5,-F(41,2400)))
OLD = (
    (3, None, 0, F(21,50), F(29,100)),
    (5, 2, 3, F(47,100), F(163,1000)),
    (7, 5, 6, F(12,25), F(23,200)),
    (9, 8, 9, F(49,100), F(89,1000)),
)
T_DEN = 10**8
NEW_RADIUS = F(47, 100)
T_LIMIT = F(1,4)
L = sum(k*abs(a) for k,a in TERMS)
K = 96
N_NEW = 448


def validate_geometry(pi):
    assert pi < ar(PI_UPPER)
    for n, _, _, r, l in OLD:
        a=pi/(2*n)
        assert 1-1/(2*a.cos())>ar(r)
        assert a.tan()/2<ar(l)
    assert F(1,2)-sum(abs(a) for _,a in TERMS)>NEW_RADIUS
    amplitude=sum((k*k-1)*abs(a) for k,a in TERMS)
    rate=sum(k*(k*k-1)*abs(a) for k,a in TERMS)
    assert F(1,2)-amplitude==F(1,100)
    assert (F(1,2)+amplitude)**2+rate**2<F(5,2)**2


def disk_points(points, n, pi, radius):
    for j in range(n):
        theta=2*pi*j/n
        append_point(points, ar(radius)*theta.cos(), ar(radius)*theta.sin())


def sqrt_upper(v):
    scaled=v*S*S
    n=isqrt(scaled.numerator//scaled.denominator)
    return F(n + (F(n*n)<scaled),S)


def old_common_polygon(center, half, pi, k=K):
    points=[]
    disk_points(points,8*k,pi,F(1,2)-F(3,S))
    assert all(x*x+y*y<=(S//2)**2 for x,y in points)
    shrink=[]
    for n, j, t, r, lip in OLD:
        rotation=center[j] if j is not None else F(0)
        angle_half=half[j] if j is not None else F(0)
        hx,hy=half[t:t+2]
        err=sqrt_upper(hx*hx+hy*hy)+lip*angle_half+F(4,S)
        factor=1-err/r
        assert 0<factor<1
        material=[]
        arc_body(material,k,pi,n,ar(rotation),arb(0),arb(0))
        tx,ty=(q*S for q in center[t:t+2])
        points.extend((round(factor*x+tx),round(factor*y+ty)) for x,y in material)
        shrink.append({"body":n,"alpha":str(factor),"support_error":str(err)})
    return hull(points),shrink


def placement_area(center, pi, new=None, k=1024):
    points=[]
    disk_points(points,8*k,pi,F(1,2))
    for n,j,t,_,_ in OLD:
        arc_body(points,k,pi,n,ar(center[j]) if j is not None else arb(0),*map(ar,center[t:t+2]))
    if new is not None:
        fourier_body(points,8*k,pi,TERMS,*map(ar,new))
    return arithmetic(hull(points),k)


def marked_hull(points):
    p=[]
    for q in sorted(points):
        if not p or p[-1][:2]!=q[:2]:
            p.append(q)
    def half(seq):
        out=[]
        for q in seq:
            while len(out)>1 and cross(out[-2],out[-1],q)<=0:
                out.pop()
            out.append(q)
        return out[:-1]
    H=half(p)+half(reversed(p))
    assert len(H)>=3 and all(cross(a,b,c)>0
        for a,b,c in zip(H,H[1:]+H[:1],H[2:]+H[:2],strict=True))
    return H


class AngleCell:
    def __init__(self,old,material,a,b,pi):
        theta=pi*ar(a+b)
        co,si=theta.cos(),theta.sin()
        rotated=[]
        for x,y in material:
            x,y=ar(F(x,S)),ar(F(y,S))
            append_point(rotated,x*co-y*si,x*si+y*co)
        error=F(4,S)+L*PI_UPPER*(b-a)
        factor=NEW_RADIUS/(NEW_RADIUS+error)
        den=S*factor.denominator
        self.scale=lcm(S,den,T_DEN)
        self.p=[(x*(self.scale//S),y*(self.scale//S),0) for x,y in old]
        multiplier=factor.numerator*(self.scale//den)
        self.q=[(x*multiplier,y*multiplier) for x,y in hull(rotated)]

    def plane(self,site):
        tx,ty=(i*(self.scale//T_DEN) for i in site)
        H=marked_hull(self.p+[(x+tx,y+ty,1) for x,y in self.q])
        v=gx=gy=0
        for (x,y,b),(u,w,c) in zip(H,H[1:]+H[:1],strict=True):
            v+=x*w-y*u
            gx+=b*w-c*y
            gy+=c*x-b*u
        return F(v-gx*tx-gy*ty,2*self.scale**2),F(gx,2*self.scale),F(gy,2*self.scale)


def bound_planes(planes,weights):
    assert len(planes)==len(weights) and all(type(w) is int and w>0 for w in weights)
    denom=sum(weights)
    b,x,y=[sum(w*r[j] for w,r in zip(weights,planes,strict=True))/denom for j in range(3)]
    return b-T_LIMIT*(abs(x)+abs(y))


def verify(old,material,pi,target,rows,denominator):
    a=F(0)
    minimum=F(1)
    assert type(denominator) is int and denominator>0
    for stop,witnesses in rows:
        assert type(stop) is int
        b=F(stop,denominator)
        assert a<b<=1
        cell=AngleCell(old,material,a,b,pi)
        assert all(type(z) is int for witness in witnesses for z in witness)
        got=bound_planes([cell.plane((x,y)) for x,y,_ in witnesses],[w for _,_,w in witnesses])
        assert got>target
        minimum=min(minimum,got)
        a=b
    assert a==1
    return minimum


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree",type=Path,help="optional original cert11d_08366 directory")
    parser.add_argument("--check",action="store_true",help="compare the archived exact results")
    args=parser.parse_args()
    ctx.prec=160
    pi=arb.pi()
    validate_geometry(pi)
    here=Path(__file__).resolve().parent
    raw=(here/"five_harmonic_inputs.json").read_bytes()
    data=json.loads(raw)
    center=list(map(F,data["center"]))
    assert len(center)==11
    material=[]
    fourier_body(material,N_NEW,pi,TERMS,arb(0),arb(0),arb(0))
    s=ar(NEW_RADIUS+T_LIMIT)
    far_area=(pi-(1/(2*s)).acos())/4+(s*s-arb(1)/4).sqrt()/2

    proofs={}
    assert set(data["cases"])=={"fixed","ball4096","ball1536"}
    for name,case in data["cases"].items():
        radius=F(case["halfwidth"])
        threshold=F(case["threshold"])
        assert radius>=0 and far_area>ar(threshold)
        old,shrink=old_common_polygon(center,[radius]*11,pi)
        rows=case["witnesses"]
        lower=verify(old,material,pi,threshold,rows,data["turn_denominator"])
        assert far_area>ar(lower)
        area=F(sum(x*v-y*u for (x,y),(u,v) in zip(old,old[1:]+old[:1],strict=True)),2*S**2)
        proofs[name]={"halfwidth":str(radius),"threshold":str(threshold),"angle_leaves":len(rows),
                      "lower":str(lower),"display":float(lower),"base_common_inner_area":str(area),
                      "base_hull_vertices":len(old),"old_shrink":shrink}

    original=placement_area(center,pi,k=2048)
    center_added=placement_area(center,pi,new=data["center_new_pose"],k=2048)
    free_pose=data["free_six_pose"]
    assert len(free_pose)==14
    free=placement_area(free_pose[:11],pi,new=free_pose[11:],k=2048)
    fixed_increment=F(proofs["fixed"]["lower"])-F(original["upper"])
    box_excluded_gap=F(proofs["ball4096"]["lower"])-F(free["upper"])
    assert F(original["upper"])<F("0.837295153")
    assert F(center_added["upper"])<F("0.838851260")
    assert F(free["upper"])<F("0.837677741")
    assert fixed_increment>F("0.0015048667") and box_excluded_gap>F("0.0000236068")

    report={
        "dependencies_sha256":{
            "five_harmonic_inputs.json":sha256(raw).hexdigest(),
            "certify_five_harmonic.py":sha256(Path(__file__).read_bytes()).hexdigest(),
            "certify_nine.py":sha256(Path(certify_nine.__file__).read_bytes()).hexdigest(),
        },
        "base_source":{"commit":data["frozen_upstream_commit"],"opt11_point_sha256":data["opt11_point_sha256"]},
        "test_harmonics":[[k,str(a)] for k,a in TERMS],
        "precision_bits":ctx.prec,
        "exact_area_results":{
            "old_material_resolution":K,"new_material_count":N_NEW,
            "cases":proofs,"fixed_original_five":original,"fixed_six_feasible":center_added,
            "free_six_feasible":free,"fixed_increase_lower":str(fixed_increment),
            "small_box_gap_above_free_six_feasible":str(box_excluded_gap),
        },
        "runtime":{"python":platform.python_version(),"python_flint":flint.__version__},
    }
    if args.tree:
        import audit_five_harmonic_tree
        report["dependencies_sha256"]["audit_five_harmonic_tree.py"]=sha256(
            Path(audit_five_harmonic_tree.__file__).read_bytes()).hexdigest()
        report["exact_tree_results"]={
            name:audit_five_harmonic_tree.audit(args.tree,center,F(data["cases"][name]["halfwidth"]))
            for name in ("ball4096","ball1536")
        }
    if args.check:
        archived=json.loads((here/"five_harmonic_results.json").read_bytes())
        assert all(v==archived["dependencies_sha256"][k] for k,v in report["dependencies_sha256"].items())
        assert all(archived[k]==v for k,v in report.items() if k not in ("dependencies_sha256","runtime"))
    print(json.dumps(report,indent=2))


if __name__=="__main__":
    main()
