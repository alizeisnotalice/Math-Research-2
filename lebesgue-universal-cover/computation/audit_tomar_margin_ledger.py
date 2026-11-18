"""Audit archived 11D lower margins, geometry outside the root and its worst recorded leaf.

Arguments: frozen v9 checkout, extracted cert11d_08366, extracted fleet11/results.
The first and later pass records are historical, except the identified worst leaf.
"""
import argparse
from fractions import Fraction
import json
from pathlib import Path
import sys

import numpy as np
from flint import arb, ctx

parser = argparse.ArgumentParser(description=__doc__)
for name in ("checkout", "certdir", "results"):
    parser.add_argument(name, type=Path)
args = parser.parse_args()
source = args.checkout.resolve()
frozen = source/"independent/full11"
sys.path.insert(0, str(frozen/"v11_frozen_8b878f23"))
import treedec11
import v11_leafcheck as verifier
ctx.prec = 128

threshold = Fraction(3,10**9)
tau = arb("0.8366")
parts = sorted(args.results.glob("*/part*.json"))
assert len(parts) == 2048
all_keys = set()
fail_keys = set()
min_saved, min_key, total, finite = 1., None, 0, 0
reference_header = None
nodes_total = 0
for path in parts:
    rec = json.loads(path.read_bytes())
    p = rec["part"]
    assert p not in all_keys
    all_keys.add(p)
    header = json.loads((args.certdir/f"part{p:04}.json").read_bytes())
    root = {k:header[k] for k in ("ns","tau","split0","lo","hi","nparts")}
    if reference_header is None:
        reference_header = root
    assert root == reference_header and header["part"] == p
    codes = np.load(args.certdir/f"part{p:04}.npy",allow_pickle=False)
    assert codes.ndim==1 and ((codes<11) | (codes>=253)).all() and codes.dtype==np.uint8
    leaf_pos = np.flatnonzero(codes>=253)
    # Number of still-unread child slots; nonempty proper prefixes, empty at the end.
    slots = 1 + np.cumsum(np.where(codes<11,1,-1),dtype=np.int64)
    assert (slots[:-1]>0).all() and slots[-1]==0
    nodes_total += len(codes)
    assert not rec["decode_errs"] and rec["leaves"]==len(leaf_pos)
    assert rec["tau"]=="0.8366"
    margins = np.load(path.with_name(path.stem+"_margins.npy"),allow_pickle=False)
    valid = np.isfinite(margins)
    assert margins.shape==(len(leaf_pos),) and margins.dtype==np.float64 and not np.isinf(margins).any()
    expected={(p,int(pos),int(codes[pos])) for pos in leaf_pos[~valid]}
    assert expected == {(p,x["pos"],x["code"]) for x in rec["fails"]}
    assert len(expected)==rec["nfail"] and valid.sum()==rec["checked"]
    fail_keys.update(expected)
    idx = int(np.nanargmin(margins))
    val = float(margins[idx])
    # Margins are nearest binary64 conversions of lower endpoints; use a full predecessor.
    assert Fraction(float(np.nextafter(val,-np.inf))) > threshold
    assert val==rec["min_margin"] and rec["min_leaf"]=={"pos":int(leaf_pos[idx]),"code":int(codes[leaf_pos[idx]])}
    if val < min_saved:
        min_saved, min_key = val,(p,int(leaf_pos[idx]),int(codes[leaf_pos[idx]]))
    total += len(margins)
    finite += int(valid.sum())

assert all_keys==set(range(2048)) and reference_header["split0"]==[2]*11
assert reference_header["ns"]==[3,5,7,9] and reference_header["nparts"]==2048
later = {}
for fname in sorted((frozen/"pass23/results").glob("*.jsonl")):
    stem=fname.stem
    suffix=stem.split("_",1)[1]
    with np.load(frozen/"pass23"/("fails_"+suffix+".npz"),allow_pickle=False) as grp:
        for line in fname.read_text().splitlines():
            item=json.loads(line)
            i=item["idx"]
            key=tuple(int(grp[k][i]) for k in ("part","pos","code"))
            assert key[:2]==(item["part"],item["pos"])
            if item["ok"] and key in fail_keys:
                margin=float(item["margin"])
                assert np.isfinite(margin) and Fraction(float(np.nextafter(margin,-np.inf)))>threshold
                later[key]=max(later.get(key,0.),margin)
assert set(later)==fail_keys and total==finite+len(fail_keys)==16193402

# Monotone area lower bound of the convex hull of the disk and a far vertex.
root_differences=[]
for i,n in enumerate(reference_header["ns"]):
    j=3*i
    w=arb(min(-reference_header["lo"][j],-reference_header["lo"][j+1],
              reference_header["hi"][j],reference_header["hi"][j+1]))
    R=1/(2*(arb.pi()/(2*n)).cos())
    s=(R**2+w**2+2*R*w*(arb.pi()/n).cos()).sqrt()
    f=(arb.pi()-(1/(2*s)).acos())/4+(s**2-arb(1)/4).sqrt()/2
    bound=f-tau
    assert bound>arb("0.000000003")
    root_differences.append(str(bound))
for r,n in [(2,5),(5,7),(8,9)]:
    assert reference_header["lo"][r]<=0 and arb(reference_header["hi"][r])>(1 if r==2 else 2)*arb.pi()/n

# Fresh 128-bit replay of the entire leaf giving the original smallest finite margin.
p,pos,code=min_key
summary,leaves=treedec11.decode_part(str(args.certdir),p,{"pos":{pos}})
assert not summary["errs"] and len(leaves)==1
tag,where,lo,hi,_,side=leaves[0]
assert (where,tag)==(pos,code) and tag==253
stats={}
ok,v,msg=verifier.check_leaf(tag,lo,hi,tau,np.array(side[0]),np.array(side[1]),0,stats)
assert ok and v-tau>arb("0.000000003")
print(json.dumps({"frozen_commit":"10af98606c98ca9ffd0d9f402cae1e2efe74b792",
    "parts":len(parts),"nodes_with_valid_prefix_structure":nodes_total,
    "recorded_leaves":total,"pass1_finite":finite,"later_resolutions":len(later),
    "recorded_pass1_min":min_saved,"recorded_pass1_min_key":min_key,
    "recorded_later_min":min(later.values()),"root_outside_margins_arb128":root_differences,
    "critical_leaf_fresh_margin_arb128":str(v-tau),"fresh_leaf_details":msg,
    "conditional_lower_rational":"836600003/1000000000",
    "scope":"All historical margin arrays, all tree code stream prefix structures, all later resolutions. Fresh arithmetic only for the root exterior and one worst archived leaf."},indent=2))
