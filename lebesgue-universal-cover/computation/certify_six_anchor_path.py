"""Certify the low-area C1--C4 six-body pose path and its Q11-excluded tube."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

from flint import arb, ctx

import certify_high_six_host_neighborhood as neighborhood
from certify_mixed_support_budget import Host,PI_UP,ar,cloud
from certify_six_structural_escape import corners

TARGET=F(".83753")
BASE=F(".83742")
CENTERS=16
AREA_PIECES=128
FIVE_SPEED=F(".000662")
POINT_SPEED=F(".003675")
OPTIONS=list(map(F,(".00037",".00036",".00035",".00034",".00033",".00032",
                    ".00030",".00028",".00026",".00024",".00022")))
EXPECTED={
    "free-p30-reoptimized-pose.json":"21e6fceeb7f2033b69a6c2a7ae633abb26b6afd3621e32f6e6aa57f9b8770593",
    "structural-six-escape-input.json":"fa4107fb115c05d559e801474457d7d985b25692b586088602930b4499f3d48b",
    "finite-point-leaf-input.json":"b99b3bf66adf904b8b99ba1e1c7f8e75abc1ea6d8e00728b9553c20b961ffebc",
    "multi-host-points.json":"9052ca388aabacc87116ba99a05d9ef98b49dfeed564e8cfd18efb647f93f282",
    "free-p30-placement.json":"438cf362c7d67245cef50900567c963e6e5cbe69671709899167c91d37c58e1d",
    "competing-fits-input.json":"31aaf17ffb9e33b4a50050614c37c78fdcfd9976360d561bd42eb1dd26070f96",
}

def sha(data):
    return sha256(data).hexdigest()

def pose_at(a,b,t):
    return [(1-t)*x+t*y for x,y in zip(a,b,strict=True)]

def midpoint(i,n):
    return F(2*i+1,2*n)

def motion(start,end,p30,pi):
    d=[y-x for x,y in zip(start,end,strict=True)]
    speeds=[]
    for n,j in ((3,-1),(5,2),(7,5),(9,8)):
        dx,dy=d[j+1:j+3]
        dist=ar(dx**2+dy**2).sqrt()
        if j>=0:
            dist+=ar(abs(d[j]))*(pi/(2*n)).tan()/2
        assert dist<ar(FIVE_SPEED)
        speeds.append(str(dist))
    q=p30["denominator"]
    for x,y in p30["points"]:
        assert F(x*x+y*y,q*q)<F(".533")**2
    v=ar(d[12]**2+d[13]**2).sqrt()+ar(F(".533")*abs(d[11]))
    assert v<ar(POINT_SPEED)
    return speeds,str(v)

def p30_clearance(host,placed):
    floor=arb(1)
    for a,b in zip(host.material,host.material[1:]+host.material[:1],strict=True):
        u,v=b[1]-a[1],a[0]-b[0]
        c=ar(F(u*a[0]+v*a[1],host.den))
        length=arb(u*u+v*v).sqrt()
        for x,y in placed:
            floor=floor.min((c-u*x-v*y)/length)
    return floor

def capsule_job(job):
    t,start,end,p30,q11,frozen=job
    ctx.prec=160
    pi=arb.pi()
    neighborhood.LOW=BASE
    position=pose_at(start,end,t)
    host=Host(position[:11],pi)
    placed=corners(position,p30,1,pi)[1]
    fixed=p30_clearance(host,placed)
    cost=F(FIVE_SPEED+POINT_SPEED,2*CENTERS)
    margin=fixed-ar(cost)
    assert margin>0
    half=F(FIVE_SPEED,2*CENTERS)
    if frozen is None:
        for d in OPTIONS:
            p=neighborhood.Caps(host,q11,pi,TARGET,d)
            try:
                cert={str(sign):p.make(sign) for sign in (1,-1)}
            except AssertionError:
                continue
            frozen={"center":str(t),"degree_radius":str(d),"hull_sha256":host.hull_hash,"hands":cert}
            break
        else:
            raise ValueError(f"No cap radius proved at t={t}")
    d=F(frozen["degree_radius"])
    assert frozen["center"]==str(t) and frozen["hull_sha256"]==host.hull_hash
    assert d in OPTIONS and d>half
    case=neighborhood.Caps(host,q11,pi,TARGET,d)
    out={"center":str(t),"center_radius":str(d),"interval_radius":str(d-half),
         "p30_center_slack":str(fixed),"p30_interval_slack":str(margin),
         "p30_interval_slack_display":float(margin),
         "augmented_host_radius_upper":case.enclosing,
         "hands":{str(s):case.check(s,frozen["hands"][str(s)]) for s in (1,-1)}}
    return out,frozen

def area_job(job):
    i,n,start,end=job
    ctx.prec=160
    pi=arb.pi()
    t=midpoint(i,n)
    host=Host(pose_at(start[:11],end[:11],t),pi)
    peri=neighborhood.perimeter_upper(host)
    err=F(FIVE_SPEED,2*n)
    e=err+host.eta
    inner=host.area
    lo=(inner-PI_UP*err**2)/(1+4*err)
    hi=inner+peri*e+PI_UP*e**2
    assert lo>BASE and hi<TARGET
    centre=inner+peri*host.eta+PI_UP*host.eta**2
    return dict(cell=i,center=str(t),inner_area=str(inner),center_upper=str(centre),
                whole_area_lower=str(lo),whole_area_upper=str(hi))

def data_files(comp):
    rows={}
    for name,expected in EXPECTED.items():
        raw=(comp/name).read_bytes()
        assert sha(raw)==expected
        rows[name]=json.loads(raw)
    return rows

def endpoint_checks(rows,start,end,pi):
    y=rows["free-p30-placement.json"]
    middle=pose_at(list(map(F,y[:11])),end[:11],F(1,2))
    c2=Host(middle,pi).area
    z=rows["competing-fits-input.json"]["zero_gain_five"]["params"]
    c3=Host(z,pi).area
    assert c2>TARGET and c3>TARGET
    def bound(x):
        h=Host(x[:11],pi)
        p=neighborhood.perimeter_upper(h)
        u=h.area+p*h.eta+PI_UP*h.eta**2
        return [str(h.area),str(u)]
    return {"C2_to_C4_straight_midpoint_lower":str(c2),
            "C2_to_C4_straight_midpoint_display":float(c2),
            "C3_anchor_lower":str(c3),"C3_anchor_display":float(c3),
            "C1_exact_pose_area_bracket":bound(start),"C4_exact_pose_area_bracket":bound(end)}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("computation",type=Path)
    p.add_argument("witness",type=Path)
    p.add_argument("--generate",action="store_true")
    p.add_argument("--workers",type=int,default=1)
    a=p.parse_args()
    assert 1<=a.workers<=16
    ctx.prec=160
    pi=arb.pi()
    assert pi<ar(PI_UP)
    rows=data_files(a.computation)
    start=list(map(F,rows["free-p30-reoptimized-pose.json"]))
    end_info=rows["structural-six-escape-input.json"]
    end=list(map(F,end_info["pose"]))
    assert len(start)==len(end)==14 and end_info["hand"]==1
    p30=rows["finite-point-leaf-input.json"]
    q11=rows["multi-host-points.json"]
    q11,d2=cloud(q11["points"],q11["point_denominator"])
    assert len(q11)==11
    speed,point=motion(start,end,p30,pi)
    r=ar(F(1,2))
    R=ar(F(".703"))
    radial=(pi-(r/R).acos())*r**2+r*(R**2-r**2).sqrt()
    assert radial>ar(TARGET)
    conf={"inputs":EXPECTED,"geometry":"linear interpolation of 14 coordinates, hand +1",
          "caps_centers":CENTERS,"area_cells":AREA_PIECES,"base":str(BASE),"target":str(TARGET),
          "whole_path_five_speed":str(FIVE_SPEED),"whole_path_P30_speed":str(POINT_SPEED)}
    old=[None]*CENTERS
    if not a.generate:
        check=json.loads(a.witness.read_bytes())
        assert check["config"]==conf
        old=check["centers"]
    cc=[(midpoint(i,CENTERS),start,end,p30,q11,f) for i,f in enumerate(old)]
    assert len(cc)==CENTERS
    ac=[(i,AREA_PIECES,start,end) for i in range(AREA_PIECES)]
    if a.workers==1:
        cres=list(map(capsule_job,cc))
        areas=list(map(area_job,ac))
    else:
        with ProcessPoolExecutor(max_workers=a.workers) as pool:
            cres=list(pool.map(capsule_job,cc))
            areas=list(pool.map(area_job,ac,chunksize=4))
    if a.generate:
        raw=json.dumps({"config":conf,"centers":[r[1] for r in cres]},separators=(",",":"))+"\n"
        temp=a.witness.with_suffix(".pending")
        temp.write_text(raw)
        temp.replace(a.witness)
    tubes=[r[0] for r in cres]
    rad=min(F(r["interval_radius"]) for r in tubes)
    lo=min(F(r["whole_area_lower"]) for r in areas)
    hi=max(F(r["whole_area_upper"]) for r in areas)
    report={"complete":True,"config":conf,"witness_sha256":sha(a.witness.read_bytes()),
            "Q11_squared_distance_numerator":d2,"analytic_body_speeds":speed,
            "analytic_p30_speed":point,"lowest_tube_radius":str(rad),
            "lowest_tube_radius_display":float(rad),"uniform_area_bracket":[str(lo),str(hi)],
            "uniform_area_display":[float(lo),float(hi)],
            "p30_min_closed_cell_clearance":min(r["p30_interval_slack_display"] for r in tubes),
            "other_references":endpoint_checks(rows,start,end,pi),
            "tube_cells":tubes,"area_cells":areas}
    print(json.dumps(report,indent=2))

if __name__=="__main__":
    main()
