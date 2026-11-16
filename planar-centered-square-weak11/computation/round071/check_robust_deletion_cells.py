#!/usr/bin/env python3
"""Independently verify R71 all-radius residual exclusion on five closed cells.

Uses only actual positive atom data and a declared rational cell configuration.
An optional production claim is compared only after independent computation.
"""
from fractions import Fraction as F
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import platform
import sys
import time


def exact(x):
    if isinstance(x, bool) or not isinstance(x, (str, int)):
        raise ValueError("only exact rational strings or integers permitted")
    return F(x)


def dist(x,p):
    return max(abs(x[0]-p[0]),abs(x[1]-p[1]))


def box_at(x,e):
    return [x[0]-e,x[0]+e,x[1]-e,x[1]+e]


def union_area(boxes):
    xs=sorted({x for b in boxes for x in b[:2]})
    area=F(0)
    for a,b in zip(xs,xs[1:]):
        mid=(a+b)/2
        ys=sorted((c,d) for l,r,c,d in boxes if l<mid<r)
        length=F(0)
        end=None
        for c,d in ys:
            if end is None or c>end:
                length+=d-c
                end=d
            elif d>end:
                length+=d-end
                end=d
        area+=(b-a)*length
    return area


def subset_cell_check(residual,box):
    """A second, non-prefix route: all 2^8-1 subsets, exact cell distance."""
    u,v,s,t=box
    minimum=None
    worst=None
    count=0
    for mask in range(1,1<<len(residual)):
        members=[residual[k] for k in range(len(residual)) if mask>>k&1]
        a=min(p[0] for p,w,i in members);b=max(p[0] for p,w,i in members)
        c=min(p[1] for p,w,i in members);d=max(p[1] for p,w,i in members)
        delta=max((b-a)/2,(d-c)/2,b-v,u-a,d-t,s-c,F(0))
        mass=sum((w for p,w,i in members),F(0))
        slack=4*delta*delta-mass
        count+=1
        if minimum is None or slack<minimum:
            minimum=slack
            worst={"subset_indices":[i for p,w,i in members],"delta":str(delta),
                   "mass":str(mass),"slack":str(slack)}
    return {"subset_count":count,"minimum_slack":str(minimum),"worst_subset":worst,
            "all_subsets_passed":minimum is None or minimum>=0}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input",type=Path,required=True)
    ap.add_argument("--config",type=Path,required=True)
    ap.add_argument("--claim",type=Path)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    start=time.perf_counter()
    raw=args.input.read_bytes(); cfgraw=args.config.read_bytes()
    data=json.loads(raw);cfg=json.loads(cfgraw)
    assert exact(data["threshold"])==1
    atoms=[(tuple(map(exact,p)),exact(w),i)for i,(p,w)in enumerate(zip(data["positions"],data["masses"]))]
    assert len(atoms)==len(data["positions"])==len(data["masses"])
    assert all(w>0 and len(p)==2 for p,w,i in atoms)
    assert len({p for p,w,i in atoms})==len(atoms)
    sint=set(data["S_indices"])
    assert sint.issubset(set(range(len(atoms))))
    residual=[a for a in atoms if a[2]not in sint]
    e=exact(cfg["epsilon"]);r0=exact(cfg["base_closed_witness_radius"])
    R=exact(cfg["common_open_witness_radius"])
    assert e>0 and R>r0+e
    centres=[tuple(map(exact,x))for x in cfg["centres"]]
    recs=[];boxes=[];certified=[]
    for k,x in enumerate(centres):
        distances={i:dist(x,p)for p,w,i in residual}
        prefixes=[]
        for d in sorted(set(distances.values())):
            inds=[i for p,w,i in residual if distances[i]<=d]
            mass=sum((w for p,w,i in residual if i in inds),F(0))
            prefixes.append({"radius":str(d),"indices":inds,"mass":str(mass),
                             "centre_slack":str(4*d*d-mass),
                             "shifted_radius":str(d-e),
                             "shifted_slack":str(4*(d-e)**2-mass)})
        # Complete prefixes include every tied atom; after the largest event,
        # mass is constant and 4r^2 increases. Before first event it is zero.
        residual_ok=all(exact(p["shifted_radius"])>0 and exact(p["shifted_slack"])>=0 for p in prefixes)
        full=[(p,w,i)for p,w,i in atoms if dist(x,p)<=r0]
        fullmass=sum((w for p,w,i in full),F(0))
        box=box_at(x,e)
        full_cell_radius=max(max(abs(z-p[j])for j in range(2)for z in box[2*j:2*j+2])for p,w,i in full)
        full_ok=full_cell_radius<R and fullmass>4*R*R
        subset_check=subset_cell_check(residual,box)
        assert residual_ok==subset_check["all_subsets_passed"]
        good=residual_ok and full_ok
        if good:
            certified.append(k);boxes.append(box)
        recs.append({"index":k,"centre":[str(v)for v in x],"prefix_count":len(prefixes),
                     "complete_closed_residual_prefixes":prefixes,
                     "residual_all_radius_exclusion_certified":residual_ok,
                     "full_witness_indices":[i for p,w,i in full],"full_witness_mass":str(fullmass),
                     "full_witness_max_radius_over_box":str(full_cell_radius),
                     "common_open_radius":str(R),"open_radius_geometric_margin":str(R-full_cell_radius),
                     "full_mass_margin_at_common_radius":str(fullmass-4*R*R),
                     "closed_box_bounds":[str(v)for v in box],"closed_box_subset_actual_deletion":good,
                     "independent_all_residual_subset_check":subset_check})
    assert certified==cfg["expected_certified_indices"]
    pairwise=[]
    for i,j in itertools.combinations(range(len(boxes)),2):
        a,b=boxes[i],boxes[j]
        separated=a[1]<b[0]or b[1]<a[0]or a[3]<b[2]or b[3]<a[2]
        assert separated
        pairwise.append({"indices":[certified[i],certified[j]],"closed_boxes_disjoint":True})
    area=union_area(boxes)
    assert area==sum(((b-a)*(d-c)for a,b,c,d in boxes),F(0))
    assert area==exact(cfg["expected_union_area"])
    claim_status=None;claimhash=None
    if args.claim:
        claimraw=args.claim.read_bytes();claimhash=hashlib.sha256(claimraw).hexdigest()
        claim=json.loads(claimraw)
        assert len(claim["receivers"])==len(recs)
        for own,other in zip(recs,claim["receivers"]):
            assert own["centre"]==other["centre"]
            a=[(p["radius"],p["mass"],p["centre_slack"])for p in own["complete_closed_residual_prefixes"]]
            b=[(p["radius"],p["mass"],p["slack"])for p in other["positive_residual_prefixes"]]
            assert a==b
            assert own["full_witness_mass"]==other["full_captured_mass"]
            if own["closed_box_subset_actual_deletion"]:
                assert other["certified_box_boundary"]=="closed"
                assert exact(other["strict_deletion_open_box_epsilon"])==e
                assert own["closed_box_bounds"]==other["strict_deletion_open_box_bounds"]
                assert exact(other["box_common_open_full_radius"])==R
                assert exact(other["box_common_open_full_margin_at_least"])==exact(own["full_mass_margin_at_common_radius"])
        assert exact(claim["strict_deletion_area_lower"])==area
        claim_status="complete prefixes, masses, five closed boxes, common open radius and union area agree"
    out={"schema":"r071-independent-robust-cell-result-v1","all_passed":True,
         "positive_atom_count":len(atoms),"residual_positive_atom_count":len(residual),
         "epsilon":str(e),"common_open_full_radius":str(R),"receivers":recs,
         "certified_indices":certified,"pairwise_disjoint_checks":pairwise,
         "certified_union_area":str(area),"certified_union_area_decimal":float(area),
         "scope":"five closed cells subset actual deletion H(mu)-H(nu); no full halo area, no P-only claim, no general C2 bound",
         "production_claim_comparison":claim_status,
         "receipt":{"python":sys.version,"platform":platform.platform(),
                    "input_sha256":hashlib.sha256(raw).hexdigest(),
                    "config_sha256":hashlib.sha256(cfgraw).hexdigest(),
                    "claim_sha256":claimhash,
                    "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    "wall_seconds":time.perf_counter()-start}}
    args.output.write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps({"all_passed":True,"certified_indices":certified,"union_area":str(area),
                      "prefix_counts":[r["prefix_count"]for r in recs],
                      "subset_checks":sum(r["independent_all_residual_subset_check"]["subset_count"]for r in recs)}))


if __name__=="__main__":
    main()
