"""Certify a real low-area escape from both published six-body Q11 tubes."""
from argparse import ArgumentParser
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json

from flint import arb, ctx

import certify_raw_path_basin as basin
from certify_high_six_host_neighborhood import perimeter_upper
from certify_mixed_support_budget import Host, ar
from certify_nine import PI_UPPER
from certify_six_anchor_path import p30_clearance, pose_at
from certify_six_structural_escape import area, corners, width_grid

WIDTH_SCALE=10**12
TARGET=F(".83753")
UPPER_GRID=1024
PATH_CELLS=64
RAW_EXPECTED=("free-p30-reoptimized-pose.json","structural-six-escape-input.json")
RAW_SHA=("21e6fceeb7f2033b69a6c2a7ae633abb26b6afd3621e32f6e6aa57f9b8770593",
         "fa4107fb115c05d559e801474457d7d985b25692b586088602930b4499f3d48b")
LINK_SPEED=F(".0011")
LINK_P30_SPEED=F(".0012")


def width_on_direction(source,angle):
    y,x=angle.sin_cos()
    half=ar(F(1,2))
    def support(nx,ny):
        trial=half
        for (tx,ty),vertices in source:
            shift=tx*nx+ty*ny
            for px,py in vertices:
                z=px*nx+py*ny
                trial=trial.max(z+shift)
                if z < -half:
                    trial=trial.max(1+z+shift)
        return trial
    return support(x,y)+support(-x,-y)


def rounded_lower(x,den):
    scaled=x*den
    a=F(int(scaled.lower().floor().unique_fmpz()),den)
    assert x>=ar(a)
    return a


def paths(comp,data):
    path1_start,path1_end=[basin.pinned(comp,{n:h})[n] for n,h in zip(RAW_EXPECTED,RAW_SHA,strict=True)]
    yield (list(map(F,path1_start)),list(map(F,path1_end["pose"])),data["six-anchor-path-audit.json"])
    yield (list(map(F,data["free-p30-placement.json"])),
           list(map(F,data["six-aligned-path-audit.json"]["rational_endpoint"])),
           data["six-aligned-path-audit.json"])


def path_diameter(p,q,audit,direction,pi):
    bounds=[]
    config=audit["config"]
    speed=F(config["five_speed"] if "five_speed" in config else config["whole_path_five_speed"])
    for j in range(PATH_CELLS):
        t=F(2*j+1,2*PATH_CELLS)
        body=corners(pose_at(p[:11],q[:11],t),None,1,pi)[0]
        bound=rounded_lower(width_on_direction(body,ar(direction)*pi),WIDTH_SCALE)
        bounds.append(bound)
    result=min(bounds)-speed/PATH_CELLS
    if "tail" in audit:
        body=corners(q[:11],None,1,pi)[0]
        endpoint=rounded_lower(width_on_direction(body,ar(direction)*pi),WIDTH_SCALE)
        result=min(result,endpoint-2*F(audit["tail"]["maximum_Hausdorff_excursion"]))
    return result,bounds


def straight_piece(job):
    start,end,p30,left,right=job
    ctx.prec=160
    pi=arb.pi()
    center=pose_at(start,end,(left+right)/2)
    host=Host(center[:11],pi)
    points=corners(center,p30,1,pi)[1]
    margin=p30_clearance(host,points)
    half=(right-left)/2
    loss=LINK_SPEED*half
    outer=host.eta+loss
    upper=host.area+perimeter_upper(host)*outer+PI_UPPER*outer**2
    accepted=upper<TARGET and margin>ar(half*(LINK_SPEED+LINK_P30_SPEED))
    return accepted,upper,str(margin-ar(half*(LINK_SPEED+LINK_P30_SPEED)))


def connection(data,end,p30,workers,pi):
    start=list(map(F,data["free-p30-placement.json"]))
    for n,j in ((3,-1),(5,2),(7,5),(9,8)):
        shift=sum((end[k]-start[k])**2 for k in (j+1,j+2))
        angular=ar(abs(end[j]-start[j]))*(pi/(2*n)).tan()/2 if j>=0 else 0
        assert ar(shift).sqrt()+angular<ar(LINK_SPEED)
    den=p30["denominator"]
    assert all(F(x*x+y*y,den*den)<F(".533")**2 for x,y in p30["points"])
    assert ar(sum((end[k]-start[k])**2 for k in (12,13))).sqrt()+ar(F(".533")*abs(end[11]-start[11]))<ar(LINK_P30_SPEED)

    # Enclose both area and fixed-hand P30 along closed parameter subintervals.
    count=64
    pending=[(F(j,count),F(j+1,count)) for j in range(count)]
    finished=[]
    with ProcessPoolExecutor(max_workers=workers) as executor:
        while pending:
            jobs=[(start,end,p30,a,b) for a,b in pending]
            nxt=[]
            for (l,r),(ok,upper,margin) in zip(pending,executor.map(straight_piece,jobs,chunksize=2),strict=True):
                if ok:
                    finished.append((l,r,upper,margin))
                else:
                    assert r-l>F(1,1024)
                    m=(l+r)/2
                    nxt.extend(((l,m),(m,r)))
            pending=nxt
    finished.sort()
    assert finished[0][0]==0 and finished[-1][1]==1
    assert all(a[1]==b[0] for a,b in zip(finished,finished[1:]))
    for l,r,upper,margin in finished:
        assert upper<TARGET and arb(margin)>0
    return {"body_speed":str(LINK_SPEED),"P30_speed":str(LINK_P30_SPEED),
            "maximum_area_upper":str(max(r[2] for r in finished)),
            "maximum_area_upper_display":float(max(r[2] for r in finished)),
            "rows":[list(map(str,row)) for row in finished]}


def main():
    p=ArgumentParser(description=__doc__)
    p.add_argument("computation",type=Path)
    p.add_argument("tree",type=Path)
    p.add_argument("pose",type=Path)
    p.add_argument("--connection-workers",type=int,default=0)
    args=p.parse_args()
    assert 0<=args.connection_workers<=12
    ctx.prec=160
    pi=arb.pi()
    assert pi<ar(PI_UPPER)
    data=basin.pinned(args.computation,basin.COMPUTATION_SHA)
    src=args.pose.read_bytes()
    spec=json.loads(src)
    assert spec["original_part"]==basin.PART and spec["original_preorder_position"]==basin.POSITION
    pose=list(map(F,spec["pose"]))
    assert spec["orientation"]==1 and len(pose)==14
    lo,hi,_=basin.raw_leaf(args.tree,pose[:11])
    assert all(a<x*basin.Q<b for a,x,b in zip(lo,pose[:11],hi,strict=True))
    p30=data["finite-point-leaf-input.json"]

    bodies,points=corners(pose,p30,1,pi)
    host=Host(pose[:11],pi)
    clearance=p30_clearance(host,points)
    assert clearance>ar(F(".00074"))
    measures,_,_=area(pose,points,pi)
    assert F(".83752367")<measures[0]<measures[1]<F(".83752378691")<TARGET

    candidate_width=F(int(width_grid((bodies,[]),UPPER_GRID,pi,False).max()),WIDTH_SCALE)
    candidate_diameter=F("1.110176")
    assert ar(candidate_width)<ar(candidate_diameter)*(pi/(2*UPPER_GRID)).cos()
    records=[]
    indices=spec["path_direction_indices"]
    assert len(indices)==2 and spec["direction_index_denominator"]==2048
    max_radius=F(0)
    for (start,end,audit),k in zip(paths(args.computation,data),indices,strict=True):
        assert isinstance(k,int) and 0<=k<2048
        assert F(audit["config"]["target"])==TARGET
        path_radius=max(F(row["center_radius"]) for row in audit["tube_cells"])
        if "tail" in audit:
            path_radius=max(path_radius,F(audit["tail"]["radius"]))
        max_radius=max(path_radius,max_radius)
        d,rows=path_diameter(start,end,audit,F(k,2048),pi)
        assert d>F("1.11126")
        records.append({"direction_in_units_pi":str(F(k,2048)),
                        "diameter_lower":str(d),"diameter_lower_display":float(d),
                        "center_width_lower":list(map(str,rows))})
    distance=(min(F(row["diameter_lower"]) for row in records)-candidate_diameter)/2
    assert max_radius==F(".00037") and distance>F(".00054")
    output={
        "complete":True,"pose_sha256":sha256(src).hexdigest(),
        "source_original_leaf":[basin.PART,basin.POSITION],"tree_sha256":basin.TREE_SHA,
        "path_source_sha256":basin.COMPUTATION_SHA|dict(zip(RAW_EXPECTED,RAW_SHA,strict=True)),
        "area_interval":list(map(str,measures)),"area_display":list(map(float,measures)),
        "fixed_P30_support_clearance":str(clearance),
        "candidate_width_upper_grid":UPPER_GRID,
        "candidate_diameter_upper":str(candidate_diameter),"both_paths":records,
        "uniform_Hausdorff_lower_any_subset_and_isometry":str(distance),
        "uniform_Hausdorff_lower_display":float(distance),
        "published_maximum_radius":str(max_radius),
        "uniform_excess_over_published_tubes":str(distance-max_radius)
    }
    if args.connection_workers:
        output["low_area_connection_to_C2"]=connection(data,pose,p30,args.connection_workers,pi)
    print(json.dumps(output,indent=2))


if __name__=="__main__":
    main()
