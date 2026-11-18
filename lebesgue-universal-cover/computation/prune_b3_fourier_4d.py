#!/usr/bin/env python3
"""Exact 4D pre-pruning for (t_B3x,t_B3y,t_Fx,t_Fy).

Only D, an inner triangle of B3, and an orientation-independent inner
octagon of the Fourier body are used.  A pruned 4D cell is therefore excluded
for every orientation of the Fourier body and every value of the remaining
16 F8 parameters.  Retained cells make no feasibility claim.
"""
from __future__ import annotations

import argparse, hashlib, itertools, json, platform, sys, time
from collections import Counter, deque
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = ROOT / "lebesgue-universal-cover/computation/outer_domain_certificate.json"
DEFAULT_B3 = ROOT / "lebesgue-universal-cover/computation/027-b3-pruning.json"
TAU = F(419, 500)

def cross(a,b,c): return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
def det(a,b): return a[0]*b[1]-a[1]*b[0]

def hull(points):
    labels={}
    for x,y,z in points: labels.setdefault((x,y),z)
    p=[(x,y,labels[(x,y)]) for x,y in sorted(labels)]
    def part(seq):
        out=[]
        for q in seq:
            while len(out)>=2 and cross(out[-2],out[-1],q)<=0: out.pop()
            out.append(q)
        return out[:-1]
    return part(p)+part(reversed(p))

def b3_inner():
    q=[(F(11,20),F(0)),(F(-11,40),F(47,100)),(F(-11,40),F(-47,100))]
    lo,hi=F(1732050807568877,10**15),F(1732050807568878,10**15)
    assert lo*lo<3<hi*hi
    ilo,ihi=1/hi,1/lo
    centers=[((ilo,ihi),(F(0),F(0))),((-ihi/2,-ilo/2),(F(1,2),F(1,2))),
             ((-ihi/2,-ilo/2),(F(-1,2),F(-1,2)))]
    upper=[]
    for x,y in q:
        row=[]
        for (xl,xh),(yl,yh) in centers:
            dx=max(abs(x-xl),abs(x-xh)); dy=max(abs(y-yl),abs(y-yh)); z=dx*dx+dy*dy
            assert z<1; row.append(str(z))
        upper.append(row)
    return q,{"sqrt3_interval":[str(lo),str(hi)],"squared_distance_upper_bounds":upper}

def fourier_inner():
    r=F(283,600)
    q=[(F(47,100),F(0)),(F(33,100),F(33,100)),(F(0),F(47,100)),
       (F(-33,100),F(33,100)),(F(-47,100),F(0)),(F(-33,100),F(-33,100)),
       (F(0),F(-47,100)),(F(33,100),F(-33,100))]
    checks=[]
    for x,y in q:
        z=x*x+y*y; assert z<r*r; checks.append(str(z))
    # h(theta)>=1/2-|1/80|-|19/1200|=283/600, hence rB is contained
    # in the support-function body, invariant under its arbitrary rotation.
    assert F(1,2)-F(1,80)-F(19,1200)==r
    return q,{"radius":str(r),"point_norm_squared":checks,
              "support_argument":"h>=1/2-1/80-19/1200=283/600; thus (283/600)B is contained for every rotation"}

def cell_record(path, box, status, reason, amin=None, hlen=None, relation=None):
    d={"path":path,"box_scaled":[[str(a),str(b)] for a,b in box],"status":status,"reason":reason}
    if amin is not None: d["minimum_twice_fan_area_scaled"]=str(amin)
    if hlen is not None: d["fan_vertices_at_center"]=hlen
    if relation is not None: d["b3_projection_relation"]=relation
    return d

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--domain",type=Path,default=DEFAULT_DOMAIN)
    ap.add_argument("--b3-certificate",type=Path,default=DEFAULT_B3)
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--max-depth",type=int,default=16)
    ap.add_argument("--node-budget",type=int,default=8000)
    a=ap.parse_args()
    assert 0<=a.max_depth<=32 and a.node_budget>0
    began=time.monotonic(); raw=a.domain.read_bytes(); inp=json.loads(raw); den=int(inp["disk_polygon_denominator"])
    # With cyclic splitting, each axis is bisected at most ceil(depth/4) times.
    assert den%(100*(2**((a.max_depth+3)//4)))==0
    disk=[(int(x),int(y),"D") for x,y in inp["disk_polygon_points"]]
    assert all(x*x+y*y<=den*den//4 for x,y,_ in disk)
    bq,bcert=b3_inner(); fq,fcert=fourier_inner()
    bq=[(int(x*den),int(y*den),f"B{i}") for i,(x,y) in enumerate(bq)]
    fq=[(int(x*den),int(y*den),f"F{i}") for i,(x,y) in enumerate(fq)]
    bbase={z:(x,y) for x,y,z in bq}; fbase={z:(x,y) for x,y,z in fq}
    # Round-14 global box: B3 coordinates +/-21/100; Fourier +/-6/25.
    root=((-21*den//100,21*den//100),(-21*den//100,21*den//100),
          (-6*den//25,6*den//25),(-6*den//25,6*den//25))
    assert all((u+v)%2==0 for u,v in root)
    threshold=(2*TAU*den*den); assert threshold.denominator==1; threshold=threshold.numerator

    # Load the prior B3-only exact leaf partition exclusively to classify
    # genuinely new coupled pruning.  It is not used in the exclusion proof.
    braw=a.b3_certificate.read_bytes(); br=json.loads(braw)
    assert int(br["coordinate_denominator"])==den
    assert tuple(map(int,br["translation_domain_scaled"]))==root[0]
    assert F(br["tau"])==TAU
    assert br["source_global_box"]["sha256"]==hashlib.sha256(raw).hexdigest()
    bn=2**int(br["max_level"])
    blo,bhi=root[0]; fine=(bhi-blo)//bn; assert fine*bn==bhi-blo
    bgrid=bytearray(bn*bn)
    for z in br["leaves"]:
        if z["status"]!="prune": continue
        x0,x1=map(int,z["x_scaled"]); y0,y1=map(int,z["y_scaled"])
        for ix in range((x0-blo)//fine,(x1-blo)//fine):
            for iy in range((y0-blo)//fine,(y1-blo)//fine): bgrid[ix*bn+iy]=1

    def reconstruct(v,t):
        x,y,tag=v
        if tag=="D": return x,y
        if tag.startswith("B"): return bbase[tag][0]+t[0],bbase[tag][1]+t[1]
        return fbase[tag][0]+t[2],fbase[tag][1]+t[3]
    def relation(box):
        (x0,x1),(y0,y1),_,_=box
        assert all((v-blo)%fine==0 for v in (x0,x1,y0,y1))
        hits=sum(bgrid[ix*bn+iy] for ix in range((x0-blo)//fine,(x1-blo)//fine)
                 for iy in range((y0-blo)//fine,(y1-blo)//fine))
        total=((x1-x0)//fine)*((y1-y0)//fine)
        return "fully_b3_projection_pruned" if hits==total else ("disjoint_from_b3_projection" if hits==0 else "mixed_with_b3_projection")
    def assess(box):
        c=tuple((u+v)//2 for u,v in box)
        h=hull(disk+[(x+c[0],y+c[1],z) for x,y,z in bq]+[(x+c[2],y+c[3],z) for x,y,z in fq])
        # A central material point must reconstruct identically.  This is the
        # explicit guard against adding centre+t twice.
        assert all(reconstruct(v,c)==(v[0],v[1]) for v in h)
        vals=[]
        for t in itertools.product(*[(u,v) for u,v in box]):
            pts=[reconstruct(v,t) for v in h]
            # Each determinant is multiaffine in t_B3x,t_B3y,t_Fx,t_Fy;
            # same-body quadratic terms cancel, and mixed B3/F terms are
            # bilinear but affine in every single coordinate.
            if any(det(p,q)<=0 for p,q in zip(pts,pts[1:]+pts[:1])):
                return "retain","fan_orientation_not_certified",None,len(h)
            vals.append(sum(det(p,q) for p,q in zip(pts,pts[1:]+pts[:1])))
        amin=min(vals)
        if amin>threshold: return "prune","strict_inner_fan_area_exceeds_tau",amin,len(h)
        return "retain","inner_fan_lower_bound_not_above_tau",amin,len(h)
    def split(box,axis):
        u,v=box[axis]; m=(u+v)//2; assert 2*m==u+v
        q=list(box); q[axis]=(u,m); left=tuple(q); q=list(box); q[axis]=(m,v); return left,tuple(q)

    queue=deque([("",root)]); leaves=[]; counts=Counter(); processed=0
    while queue and processed<a.node_budget:
        path,box=queue.popleft(); processed+=1; stat,why,amin,hlen=assess(box); counts[(stat,why)]+=1
        if stat=="prune": leaves.append(cell_record(path,box,stat,why,amin,hlen,relation(box)))
        elif len(path)<a.max_depth:
            left,right=split(box,len(path)%4); queue.extend(((path+"0",left),(path+"1",right)))
        else: leaves.append(cell_record(path,box,stat,why,amin,hlen,relation(box)))
    while queue:
        path,box=queue.popleft(); counts[("retain","node_budget_not_processed")]+=1
        leaves.append(cell_record(path,box,"retain","node_budget_not_processed",relation=relation(box)))
    pr=[z for z in leaves if z["status"]=="prune"]
    def vol(rows):
        return sum(__import__('functools').reduce(lambda x,y:x*y,(int(v)-int(u) for u,v in z["box_scaled"]),1) for z in rows)
    domainvol=__import__('functools').reduce(lambda x,y:x*y,(v-u for u,v in root),1)
    rel=Counter(z.get("b3_projection_relation","not_pruned") for z in pr)
    relvol={k:str(F(vol([z for z in pr if z.get("b3_projection_relation")==k]),domainvol)) for k in sorted(rel)}
    out={"status":"completed_exact_4d_projection_pruning",
         "scope":"Strict only for the four translations t_B3,t_F. Pruned cells are excluded uniformly in Fourier rotation and the other 16 F8 parameters; retained cells are unresolved, not feasible. This is not a 20D lower-bound proof.",
         "tau":str(TAU),"coordinate_denominator":str(den),"root_box":"t_B3 in [-21/100,21/100]^2; t_F in [-6/25,6/25]^2",
         "source_global_box":{"path":str(a.domain),"sha256":hashlib.sha256(raw).hexdigest()},
         "comparison_b3_certificate":{"path":str(a.b3_certificate),"sha256":hashlib.sha256(braw).hexdigest(),"role":"classification only, not used to prove any 4D leaf"},
         "kernel":{"disk":"512 rational points inside fixed radius-1/2 disk","B3_triangle":[[str(F(x,den)),str(F(y,den))] for x,y,_ in bq],"B3_membership":bcert,
                   "Fourier_octagon":[[str(F(x,den)),str(F(y,den))] for x,y,_ in fq],"Fourier_membership":fcert},
         "criterion":{"orientation":"all 16 4D cell corners have strict positive consecutive fan determinants","area":"twice shoelace fan area is multiaffine, so its cell minimum is a corner","center_guard":"tag-to-base reconstruction at center equals stored centered hull vertex"},
         "max_depth":a.max_depth,"node_budget":a.node_budget,"processed_nodes":processed,"budget_exhausted":processed>=a.node_budget,
         "counts_by_assessment":{f"{x}:{y}":n for (x,y),n in sorted(counts.items())},"terminal_leaf_counts":{"all":len(leaves),"pruned":len(pr),"retained":len(leaves)-len(pr)},
         "four_dimensional_volume_fraction":{"meaning":"fractions of the original 4D translation box, not an absolute R^4 volume","pruned":str(F(vol(pr),domainvol)),"retained":str(F(domainvol-vol(pr),domainvol)),
             "pruned_by_relation_to_B3_projection_counts":dict(sorted(rel.items())),"pruned_by_relation_to_B3_projection_volume_fraction":relvol},
         "leaves":leaves,"python":sys.version,"platform":platform.platform(),"script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"elapsed_seconds":time.monotonic()-began}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"output":str(a.output),"processed_nodes":processed,"budget_exhausted":out["budget_exhausted"],"terminal_leaf_counts":out["terminal_leaf_counts"],"four_dimensional_volume_fraction":out["four_dimensional_volume_fraction"],"elapsed_seconds":out["elapsed_seconds"]},indent=2,sort_keys=True))
if __name__=="__main__": main()
