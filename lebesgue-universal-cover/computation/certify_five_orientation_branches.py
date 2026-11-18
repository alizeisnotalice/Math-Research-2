"""Classify the 501 exceptional higher-target rectangles by the invariant ninth-body angle."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np
from flint import arb,ctx
from scipy.sparse.csgraph import connected_components

sys.path.insert(0,str(Path(__file__).resolve().parent))
from tomar_leaf_box import leaf_box
from verify_tomar_fan_table import PINNED_MANIFEST

def digest(path):
    return sha256(path.read_bytes()).hexdigest()

def main():
    if not __debug__:raise SystemExit("Assertions required.")
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cert",type=Path)
    parser.add_argument("high",type=Path)
    args=parser.parse_args()
    high=args.high/"five-global-localization.json"
    basin=args.high/"five-body-contact-basin.json"
    assert digest(high)=="d85941b19b47cda83fc86268daa50569db2d80118b7fc722d4bf90adfa457504"
    assert digest(basin)=="282fafa2a4d37f72871c082b9e3e3e71c19fdb9ee4889cfa1a08fef19dcd150a"
    manifest=args.cert/"SHA256SUMS"
    assert digest(manifest)==PINNED_MANIFEST
    pins={n.removeprefix("./"):h for h,n in (s.split() for s in manifest.read_text().splitlines())}
    S=2**44;ids=[];lows=[];highs=[]
    for k,indices in json.loads(high.read_bytes())["retained"].items():
        part=int(k);head=args.cert/f"part{part:04}.json";tree=head.with_suffix(".npy")
        assert digest(head)==pins[head.name] and digest(tree)==pins[tree.name]
        meta=json.loads(head.read_bytes());codes=np.load(tree,allow_pickle=False)
        prefix=np.r_[1,1+np.cumsum(np.where(codes<11,1,-1))]
        for pos in indices:
            c,w=leaf_box(prefix,codes,part,pos,meta)
            lo=[(F(x)-F(y))*S for x,y in zip(c,w)];hi=[(F(x)+F(y))*S for x,y in zip(c,w)]
            assert all(v.denominator==1 for v in lo+hi)
            ids.append([part,pos]);lows.append(list(map(int,lo)));highs.append(list(map(int,hi)))
    lo=np.array(lows,np.int64);hi=np.array(highs,np.int64)
    assert lo.shape==hi.shape==(501,11)
    adjacency=np.stack([np.all(np.maximum(x,lo)<=np.minimum(y,hi),axis=1) for x,y in zip(lo,hi)])
    assert np.all(adjacency==adjacency.T)
    n,group=connected_components(adjacency,directed=False)
    comps=sorted((np.flatnonzero(group==i).tolist()for i in range(n)),key=len,reverse=True)
    assert list(map(len,comps))==[346,130,25]
    ctx.prec=160
    period=2*arb.pi()/9
    for i in comps[0]:
        assert arb(int(lo[i,8]))/S>arb(3)/10 and period-arb(int(hi[i,8]))/S>arb(3)/10
    for i in comps[1]+comps[2]:
        assert arb(int(hi[i,8]))/S<arb(9)/200 or period-arb(int(lo[i,8]))/S<arb(9)/200
    seam=[]
    keep=[i for i in range(11)if i!=8]
    for i in comps[1]:
        if not (arb(int(lo[i,8]))<=period*S<=arb(int(hi[i,8]))):continue
        for j in comps[2]:
            if not lo[j,8]<=0<=hi[j,8]:continue
            if np.all(np.maximum(lo[i,keep],lo[j,keep])<=np.minimum(hi[i,keep],hi[j,keep])):
                seam.append([ids[i],ids[j]])
    assert len(seam)==1324
    data=json.loads(basin.read_bytes());cent=list(map(F,data["raw_center"]));width=list(map(F,data["raw_halfwidth"]))
    left=[(c-w)*S for c,w in zip(cent,width)];right=[(c+w)*S for c,w in zip(cent,width)]
    meets=[ids[i] for i in range(len(ids)) if all(int(lo[i,j])<=right[j] and left[j]<=int(hi[i,j]) for j in range(11))]
    assert sorted(meets)==[[176,35228],[176,35845],[180,23955],[180,24887]]
    print(json.dumps(dict(parent_sha256=digest(high),basin_sha256=digest(basin),input_manifest_sha256=PINNED_MANIFEST,
        rectangles=len(ids),rectangular_component_sizes=list(map(len,comps)),ninth_periodic_seam_pairs=len(seam),
        seam_example=seam[0],near_raw_basin_intersections=meets,
        angle_distance_to_ninth_lattice={"near_strict_lower":"3/10","far_strict_upper":"9/200"},
        components=[[ids[i]for i in comp] for comp in comps]),separators=(",",":")))

if __name__=="__main__":main()
