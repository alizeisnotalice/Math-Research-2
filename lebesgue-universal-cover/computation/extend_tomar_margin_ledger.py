"""Recheck all recorded weak 11D leaves at a given higher decimal threshold.

Run after audit_tomar_margin_ledger.py on the same inputs. Example target 0.836601.
"""
import argparse
from collections import defaultdict
from decimal import Decimal
from fractions import Fraction
import json
from pathlib import Path
import sys

import numpy as np
from flint import arb

parser = argparse.ArgumentParser(description=__doc__)
for name in ("checkout", "certdir", "results"):
    parser.add_argument(name,type=Path)
parser.add_argument("target",type=Decimal)
parser.add_argument("--refinement-depth",type=int,default=0)
args=parser.parse_args()
source=args.checkout.resolve()
sys.path.insert(0,str(source/"checker"))
import lbcheck
sys.path.insert(0,str(Path(__file__).resolve().parent))
from tomar_leaf_box import leaf_box

old=Decimal("0.8366")
target=str(args.target)
assert args.target>old and repr(float(target))==target
gap=Fraction(args.target-old)
later={}
for f in (source/"independent/full11/pass23/results").glob("*.jsonl"):
    for line in f.read_text().splitlines():
        x=json.loads(line)
        if x["ok"]:
            key=(x["part"],x["pos"])
            later[key]=max(later.get(key,0.),x["margin"])

weak=defaultdict(list)
total=0
for recfile in sorted(args.results.glob("*/part*.json")):
    part=json.loads(recfile.read_bytes())["part"]
    vals=np.load(recfile.with_name(recfile.stem+"_margins.npy"))
    codes=np.load(args.certdir/f"part{part:04}.npy")
    positions=np.flatnonzero(codes>=253)
    uncertain=np.flatnonzero(~(np.nextafter(vals,-np.inf)>np.nextafter(float(gap),np.inf)))
    for rank in uncertain:
        pos=int(positions[rank]); val=float(vals[rank])
        if np.isnan(val):
            val=later[part,pos]
        if Fraction(float(np.nextafter(val,-np.inf)))<=gap:
            weak[part].append((pos,int(codes[pos]),val))
    total+=len(vals)
assert total==16193402

rows=[]
for part,items in sorted(weak.items()):
    meta=json.loads((args.certdir/f"part{part:04}.json").read_bytes())
    codes=np.load(args.certdir/f"part{part:04}.npy")
    prefix=np.r_[1,1+np.cumsum(np.where(codes<11,1,-1))]
    with np.load(args.certdir/f"part{part:04}_multi.npz") as side:
        positions=side["pos"]; offsets=np.r_[0,np.cumsum(side["nstruct"])]
        shifts,lam=side["shifts"],side["lam"]
        assert np.array_equal(positions,np.flatnonzero(codes==253))
        for pos,code,old_margin in items:
            c,h=leaf_box(prefix,codes,part,pos,meta)
            ch=lbcheck.Checker(meta["ns"],float(target),meta["N"],meta["eps"])
            slices=None
            if code==253:
                i=int(np.searchsorted(positions,pos)); assert positions[i]==pos
                sl=slice(int(offsets[i]),int(offsets[i+1])); assert len(shifts[sl])==len(lam[sl])
                slices=(shifts[sl],lam[sl])
            calls=[0]
            max_depth=[0]
            def prove(cc,hh,depth):
                calls[0]+=1
                max_depth[0]=max(max_depth[0],depth)
                if code==255:
                    ok=ch.check_far_leaf(cc-hh,cc+hh,meta["T"]);msg="corner"
                elif code==254:
                    ok,msg=ch.check_bound_leaf(cc,hh)
                else:
                    ok,msg=ch.check_multi_leaf(cc,hh,*slices)
                    if not ok and depth>0:
                        # Preserve the parent's sample positions when a smaller box
                        # changes the hull combinatorics at its newly chosen centre.
                        relative=(c-cc+slices[0]*h)/hh
                        ok,msg=ch.check_multi_leaf(cc,hh,relative,slices[1])
                    if not ok and depth>0:
                        ok,msg=ch.check_bound_leaf(cc,hh)
                if ok:
                    return True,msg
                if depth>=args.refinement_depth:
                    return False,msg
                # Redraw the witnesses on both exact half-boxes. A failed fan
                # of the parent gives no constraint to either subproblem.
                k=int(np.argmax(hh*np.array([1,1,.5,1,1,.5,1,1,.5,1,1])))
                half=hh.copy();half[k]/=2
                assert Fraction(half[k])*2==Fraction(hh[k])
                for sign in (-1,1):
                    child=cc.copy();child[k]+=sign*half[k]
                    assert Fraction(child[k])==Fraction(cc[k])+sign*Fraction(half[k])
                    ok,detail=prove(child,half,depth+1)
                    if not ok:
                        return False,detail
                return True,"both dyadic children passed"
            ok,msg=prove(c,h,0)
            rows.append({"part":part,"pos":pos,"code":code,"old_margin":old_margin,
                         "passed":bool(ok),"boxes_checked":calls[0],"depth":max_depth[0],"detail":msg})

meta=json.loads((args.certdir/"part0000.json").read_bytes())
root=[]
for i,n in enumerate(meta["ns"]):
    r=1-1/(2*(arb.pi()/(2*n)).cos())
    j=3*i
    w=arb(min(-meta["lo"][j],meta["hi"][j],-meta["lo"][j+1],meta["hi"][j+1]))
    lower=arb.pi()*(arb(1)/4+r*r)/2+w*(arb(1)/2+r)
    assert lower>arb(target)
    root.append(str(lower))
print(json.dumps({"target":target,"historical_leaves":total,"weak":len(rows),"refinement_depth":args.refinement_depth,
    "fresh_passes":sum(r["passed"] for r in rows),
    "fresh_failures":sum(not r["passed"] for r in rows),
    "root_outside_semicircle_bounds":root,"rows":rows},indent=2))
