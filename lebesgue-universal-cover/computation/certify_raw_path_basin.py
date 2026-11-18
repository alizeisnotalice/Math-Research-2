"""Authenticate an original 11D leaf and certify a full-dimensional path basin within it."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json

from flint import arb, ctx
import numpy as np

from certify_six_anchor_path import p30_clearance, pose_at
from certify_mixed_support_budget import Host, PI_UP, ar
from certify_high_six_host_neighborhood import perimeter_upper
from certify_six_structural_escape import corners

PART, POSITION = 1205, 152529
Q = 2**44
TARGET = F(".83753")
T = F(1,2)
RAW_CORNER = "00110101001"
COMPUTATION_SHA = {
    "six-anchor-path-audit.json":"48c08051ddbe7104e271b1073eb33b062d36b24566bb02fd01159bdeb511ad96",
    "six-aligned-path-audit.json":"1cdcfc37825ff0c0a46da430ba3c87b33545668d2d6943a77f97645ed3543ff4",
    "free-p30-placement.json":"438cf362c7d67245cef50900567c963e6e5cbe69671709899167c91d37c58e1d",
    "finite-point-leaf-input.json":"b99b3bf66adf904b8b99ba1e1c7f8e75abc1ea6d8e00728b9553c20b961ffebc",
}
TREE_SHA = {
    "part1205.json":"0529338eac7b0fe9e02577cb9d2de7bbe597a3f9ae49dd96561a3c93bb8d1b4c",
    "part1205.npy":"df2413c532134b4853c744da6e2fecac592c106bc81d62f9832c78b0cbe96112",
}


def pinned(base, mapping):
    files={}
    for name,digest in mapping.items():
        raw=(base/name).read_bytes()
        assert sha256(raw).hexdigest()==digest
        if name.endswith(".json"):
            files[name]=json.loads(raw)
    return files


def raw_leaf(tree_root, anchor):
    data=pinned(tree_root,TREE_SHA)["part1205.json"]
    assert data["ns"]==[3,5,7,9] and data["nparts"]==2048
    assert data["split0"]==[2]*11 and data["part"]==PART
    bounds=[[F(v)*Q for v in data[k]] for k in ("lo","hi")]
    assert all(v.denominator==1 for row in bounds for v in row)
    lo,hi=[list(map(int,row)) for row in bounds]
    for j in range(11):
        assert (lo[j]+hi[j])%2==0
        mid=(lo[j]+hi[j])//2
        if (PART >> (10-j))&1:
            lo[j]=mid
        else:
            hi[j]=mid
    raw=np.load(tree_root/"part1205.npy",allow_pickle=False)
    assert raw.dtype==np.uint8 and raw.ndim==1 and len(raw)==data["boxes"]
    assert np.all((raw<11)|(raw>=253))
    slots=np.r_[1,1+np.cumsum(np.where(raw<11,1,-1))]
    assert slots[-1]==0 and np.all(slots[:-1]>0)
    k=0
    steps=[]
    while raw[k]<11:
        axis=int(raw[k])
        assert (lo[axis]+hi[axis])%2==0
        m=(lo[axis]+hi[axis])//2
        if anchor[axis]*Q<m:
            hi[axis]=m
            k+=1
            steps.append([axis,0])
        else:
            lo[axis]=m
            k=k+2+int(np.flatnonzero(slots[k+2:]==slots[k])[0])
            steps.append([axis,1])
    assert k==POSITION and raw[k]==253
    assert all(l<x*Q<h for l,x,h in zip(lo,anchor,hi,strict=True))
    return lo,hi,steps


def excursion(lo,hi,anchor,pi):
    delta=[max(x-F(a,Q),F(b,Q)-x) for a,b,x in zip(lo,hi,anchor,strict=True)]
    values=[]
    for n,j in ((3,-1),(5,2),(7,5),(9,8)):
        velocity=ar(delta[j+1]**2+delta[j+2]**2).sqrt()
        if j>=0:
            velocity+=ar(delta[j])*(pi/(2*n)).tan()/2
        values.append(velocity)
    scale=10**12
    radius=F(max(int(float(v.upper())*scale)+2 for v in values),scale)
    assert all(v<ar(radius) for v in values)
    return radius


def select_child(lo,hi,anchor,acceptable,pi):
    lo,hi=lo.copy(),hi.copy()
    trace=[]
    last=None
    while True:
        r=excursion(lo,hi,anchor,pi)
        if acceptable(r):
            assert last is not None
            return lo,hi,trace,r,last
        last=r
        size=[F(b-a,6 if j in (2,5,8) else 1)
              for j,(a,b) in enumerate(zip(lo,hi,strict=True))]
        j=max(range(11),key=size.__getitem__)
        assert (hi[j]+lo[j])%2==0 and lo[j]<hi[j]
        m=(hi[j]+lo[j])//2
        if anchor[j]*Q<m:
            hi[j]=m
            trace.append([j,0])
        else:
            lo[j]=m
            trace.append([j,1])


def main():
    ap=ArgumentParser(description=__doc__)
    ap.add_argument("computation",type=Path)
    ap.add_argument("tree",type=Path)
    args=ap.parse_args()
    ctx.prec=160
    pi=arb.pi()
    assert pi<ar(PI_UP)
    data=pinned(args.computation,COMPUTATION_SHA)
    aligned=data["six-aligned-path-audit.json"]
    first=data["six-anchor-path-audit.json"]
    assert all(r["complete"] and F(r["uniform_area_bracket"][1])<TARGET and
               F(r["config"]["target"])==TARGET for r in (first,aligned))
    start=list(map(F,data["free-p30-placement.json"]))
    end=list(map(F,aligned["rational_endpoint"]))
    anchor=pose_at(start,end,T)
    q=F(aligned["config"]["five_speed"])
    at_t=max(F(row["center_radius"])-q*abs(F(row["center"])-T)
             for row in aligned["tube_cells"] if abs(F(row["center"])-T)<=F(1,64))
    assert at_t>0
    p30=data["finite-point-leaf-input.json"]
    lo,hi,source_trace=raw_leaf(args.tree,anchor[:11])

    # An eleven-dimensional thickening of the rational straight portion stays in the leaf.
    offset=[F(".0000012" if j in (2,5,8) else ".0000006") for j in range(11)]
    assert all(F(a,Q)<min(x,y)-d<=max(x,y)+d<F(b,Q)
               for a,x,y,b,d in zip(lo,start[:11],end[:11],hi,offset,strict=True))
    transverse=F(".00000105")
    assert (pi/10).tan()/2<ar(F(1,6))
    assert ar(F(".0000006"))*arb(2).sqrt()+ar(F(".0000012")/6)<ar(transverse)
    assert all(arb(row["p30_interval_slack"])>ar(transverse) for row in aligned["tube_cells"])
    assert F(aligned["minimum_tube_radius"])>transverse
    sweep_lo,sweep_hi=map(F,aligned["uniform_area_bracket"])
    sweep_floor=(sweep_lo-PI_UP*transverse**2)/(1+4*transverse)
    sweep_ceil=sweep_hi*(1+4*transverse)+PI_UP*transverse**2
    assert F(".83742")<sweep_floor<sweep_ceil<TARGET

    # An authenticated corner gives an isometry-invariant failure of raw-box capture.
    corner=[F((b if bit=="1" else a),Q)
            for bit,a,b in zip(RAW_CORNER,lo,hi,strict=True)]
    large=Host(corner,pi).area
    max_path_area=max(F(row["uniform_area_bracket"][1]) for row in (first,aligned))
    max_tube_radius=max(F(row["center_radius"]) for audit in (first,aligned)
                        for row in audit["tube_cells"])
    assert F(aligned["tail"]["radius"])<max_tube_radius
    separation=F(".00187")
    assert large>max_path_area*(1+4*separation)+PI_UP*separation**2
    assert large>TARGET and separation>max_tube_radius

    center=Host(anchor[:11],pi)
    points=corners(anchor,p30,1,pi)[1]
    margin=p30_clearance(center,points)
    assert margin>0
    peri=perimeter_upper(center)

    def area_upper(r):
        e=center.eta+r
        return center.area+peri*e+PI_UP*e**2

    def fits_tube(r):
        return r<at_t and ar(r)<margin

    def qualifies(r):
        return area_upper(r)<TARGET and fits_tube(r)

    tlo,thi,ttrace,tr,_=select_child(lo,hi,anchor[:11],fits_tube,pi)
    blo,bhi,trace,r,last=select_child(lo,hi,anchor[:11],qualifies,pi)
    # Preserve a strict open margin around the entire closed descendant.
    open_step=F(".00000001")
    assert qualifies(r+open_step)
    lower=(center.area-PI_UP*r**2)/(1+4*r)
    assert lower>F(".83742")
    assert all(a<x*Q<b for a,b,x in zip(blo,bhi,anchor[:11],strict=True))
    print(json.dumps({
        "complete":True,"computation_sha256":COMPUTATION_SHA,"tree_sha256":TREE_SHA,
        "part":PART,"raw_position":POSITION,"coordinate_denominator":Q,
        "raw_source_steps":source_trace,"raw_leaf_endpoints":[lo,hi],
        "whole_rational_aligned_path_in_raw_leaf":True,
        "full_dimensional_connected_sweep":{"coordinate_halfwidths":list(map(str,offset)),
            "base_path_Hausdorff_excursion_upper":str(transverse),
            "area_interval":[str(sweep_floor),str(sweep_ceil)],
            "area_display":[float(sweep_floor),float(sweep_ceil)],
            "unused_path_tube_radius":str(F(aligned["minimum_tube_radius"])-transverse)},
        "path_parameter":str(T),"anchored_pose":[str(x) for x in anchor],
        "raw_obstruction":{"corner_bits":RAW_CORNER,"area_lower":str(large),
            "area_lower_display":float(large),"all_path_host_distance_lower":str(separation),
            "minimum_excess_over_published_tubes":str(separation-max_tube_radius)},
        "tube_descendant":{"steps":ttrace,"endpoints":[tlo,thi],
            "relative_parameter_volume":str(F(1,2**len(ttrace))),
            "path_center_Hausdorff_excursion_upper":str(tr),
            "unused_path_tube_radius":str(at_t-tr),
            "fixed_P30_margin":str(margin-ar(tr))},
        "connected_sublevel_descendant":{"steps":trace,"endpoints":[blo,bhi],
            "relative_parameter_volume":str(F(1,2**len(trace))),
            "area_interval":[str(lower),str(area_upper(r))],
            "area_display":[float(lower),float(area_upper(r))],
            "path_center_Hausdorff_excursion_upper":str(r),
            "parent_Hausdorff_excursion_upper":str(last),
            "unused_path_tube_radius":str(at_t-r),
            "open_Hausdorff_neighborhood":str(open_step),
            "fixed_P30_margin":str(margin-ar(r))}},indent=2))


if __name__=="__main__":
    main()
