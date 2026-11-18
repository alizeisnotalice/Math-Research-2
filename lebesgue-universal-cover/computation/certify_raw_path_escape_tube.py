"""Verify the Q11 tube and fixed rational X18 fits on the connection to the raw-leaf escape host."""
from argparse import ArgumentParser
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

from flint import arb, ctx

import certify_high_six_host_neighborhood as neighborhood
from certify_mixed_support_budget import Host, ar, cloud
from certify_six_anchor_path import pose_at
from certify_x18_lower_hosts import inner_hull, positivity

TARGET=F(".83753")
BASE=F(".83742")
MIN_RADIUS=F(".000356")
ESCAPE_RADIUS=F(".00042")
CURVED_SPEED=F(".0011")
X18_CELLS=16
X18_ARCS=128
INPUTS={
    "raw-path-escape-pose.json":"a6600ff20f1b0f1fc3b3a9edd29a260d478a86dbee3cd75c6aec0c2e2e849904",
    "raw-path-escape-connection-audit.json":"cafb4cfdcb29dd34a2a46e621721dfe7a3e02b4496a2b360e2741b9340da9de3",
    "free-p30-placement.json":"438cf362c7d67245cef50900567c963e6e5cbe69671709899167c91d37c58e1d",
    "multi-host-points.json":"9052ca388aabacc87116ba99a05d9ef98b49dfeed564e8cfd18efb647f93f282",
    "widened-projected-obstruction-18-points.json":"65e9a5b1aa78e796b0910303d04b4c16a18af403640ddf6e322b7ab7405af71a",
}


def digest(raw):
    return sha256(raw).hexdigest()


def load_inputs(comp):
    rows={}
    for name,expected in INPUTS.items():
        raw=(comp/name).read_bytes()
        assert digest(raw)==expected
        rows[name]=json.loads(raw)
    return rows


def cap_job(job):
    start,end,q11,t,row=job
    ctx.prec=160
    pi=arb.pi()
    neighborhood.LOW=BASE
    assert row["t"]==str(t)
    host=Host(pose_at(start[:11],end[:11],t),pi)
    assert host.hull_hash==row["hull_sha256"]
    radius=F(row["radius"])
    assert 0<radius<=ESCAPE_RADIUS
    proof=neighborhood.Caps(host,q11,pi,TARGET,radius)
    out={str(hand):proof.check(hand,row["hands"][str(hand)]) for hand in (1,-1)}
    return {"t":str(t),"radius":str(radius),"enclosing_disk_radius":proof.enclosing,
            "hands":out}


def x18_fit(start,end,cloud,poses,t):
    p=pose_at(start[:11],end[:11],t)
    poly,den=inner_hull(list(map(ar,p)),X18_ARCS)
    return {str(w["sign"]):positivity(poly,den,w,cloud) for w in poses}


def main():
    ap=ArgumentParser(description=__doc__)
    ap.add_argument("computation",type=Path)
    ap.add_argument("witness",type=Path)
    ap.add_argument("--workers",type=int,default=1)
    args=ap.parse_args()
    assert 1<=args.workers<=12
    ctx.prec=160
    pi=arb.pi()
    r=ar(F(1,2))
    R=ar(neighborhood.R)
    assert (pi-(r/R).acos())*r**2+r*(R**2-r**2).sqrt()>ar(TARGET)

    data=load_inputs(args.computation)
    archived=data["raw-path-escape-connection-audit.json"]
    assert archived["complete"] and archived["pose_sha256"]==INPUTS["raw-path-escape-pose.json"]
    connected=archived["low_area_connection_to_C2"]
    assert F(connected["body_speed"])==CURVED_SPEED
    assert F(connected["maximum_area_upper"])<TARGET
    link_rows=connected["rows"]
    assert len(link_rows)==291 and F(link_rows[0][0])==0 and F(link_rows[-1][1])==1
    assert all(F(a[1])==F(b[0]) for a,b in zip(link_rows,link_rows[1:]))
    assert all(F(l)<F(r) and F(upper)<TARGET and arb(margin)>0 for l,r,upper,margin in link_rows)
    start=list(map(F,data["free-p30-placement.json"]))
    end=list(map(F,data["raw-path-escape-pose.json"]["pose"]))

    qdata=data["multi-host-points.json"]
    q11,qdiam=cloud(qdata["points"],qdata["point_denominator"])
    assert len(q11)==11
    x18=data["widened-projected-obstruction-18-points.json"]
    assert len(x18["points"])==18
    pdiam=max((x-a)**2+(y-b)**2 for x,y in x18["points"] for a,b in x18["points"])
    assert pdiam==x18["diameter_squared_numerator"] and pdiam<x18["denominator"]**2

    raw=args.witness.read_bytes()
    w=json.loads(raw)
    assert F(w["uniform_radius_requested"])==MIN_RADIUS and w["inputs"]==INPUTS
    assert F(w["endpoint"]["radius"])==ESCAPE_RADIUS
    # The frozen Q11 angular covers may use finer phase pieces than the
    # independently verified 291 area pieces. Both partitions cover [0,1].
    phase=[]
    jobs=[]
    previous=F(0)
    for row in w["phase_rows"]:
        i,k=row["cell"]
        assert 0<=i<2**k and 0<=k<=14
        left,right=F(i,2**k),F(i+1,2**k)
        assert left==previous
        previous=right
        mid=(left+right)/2
        radius=F(row["radius"])-CURVED_SPEED*(right-left)/2
        assert radius>=MIN_RADIUS
        phase.append((left,right,radius))
        jobs.append((start,end,q11,mid,row))
    assert previous==1
    link_radius=min(row[2] for row in phase)
    jobs.append((start,end,q11,F(1),w["endpoint"]))

    for key in ("endpoint_poses","whole_path_poses"):
        assert [v["sign"] for v in w["x18"][key]]==[1,-1]
    endpoint_x18=x18_fit(start,end,x18,w["x18"]["endpoint_poses"],F(1))
    assert min(endpoint_x18.values())>ESCAPE_RADIUS
    x18_margins={str(sign):[] for sign in (1,-1)}
    for i in range(X18_CELLS):
        t=F(2*i+1,2*X18_CELLS)
        for sign,value in x18_fit(start,end,x18,w["x18"]["whole_path_poses"],t).items():
            x18_margins[sign].append(value-CURVED_SPEED/(2*X18_CELLS))
    x18_clearance={sign:min(row) for sign,row in x18_margins.items()}
    assert min(x18_clearance.values())>max(F(row["radius"]) for row in w["phase_rows"])

    if args.workers==1:
        certificates=list(map(cap_job,jobs))
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as workers:
            certificates=list(workers.map(cap_job,jobs))
    print(json.dumps({
        "complete":True,"witness_sha256":digest(raw),"inputs":INPUTS,
        "target_area":str(TARGET),"base_area":str(BASE),
        "Q11_max_squared_distance_numerator":qdiam,"X18_max_squared_distance_numerator":pdiam,
        "link_source_area_cells":len(link_rows),"link_area_upper":connected["maximum_area_upper"],
        "escape_area_interval":archived["area_interval"],"escape_endpoint":certificates[-1],
        "path_uniform_Q11_tube_radius":str(link_radius),
        "path_phase_cells":[dict(left=str(l),right=str(r),interval_radius=str(d),**row)
                            for (l,r,d),row in zip(phase,certificates[:-1],strict=True)],
        "X18_endpoint_clearance":{sign:str(value) for sign,value in endpoint_x18.items()},
        "X18_endpoint_clearance_display":{sign:float(value) for sign,value in endpoint_x18.items()},
        "X18_path_support_clearance":{sign:str(value) for sign,value in x18_clearance.items()},
        "X18_path_support_clearance_display":{sign:float(value) for sign,value in x18_clearance.items()},
        "X18_path_phase_margins":{sign:list(map(str,row)) for sign,row in x18_margins.items()},
        "X18_path_phase_count":X18_CELLS,
    },indent=2))


if __name__=="__main__":
    main()
