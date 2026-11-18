"""Bind the later-pass inputs to the original 11D tree slots, boxes and side data."""
from collections import defaultdict
from fractions import Fraction
import json
from pathlib import Path
import sys

import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
from tomar_leaf_box import leaf_box

check,cert,results=map(Path,sys.argv[1:])
later=defaultdict(list)
for fn in sorted((check/"independent/full11/pass23").glob("fails_*.npz")):
    with np.load(fn,allow_pickle=False) as snap:
        d={k:snap[k] for k in snap.files}
        assert str(d["tau_str"])=="0.8366"
        for j,(p,i,code) in enumerate(zip(d["part"],d["pos"],d["code"])):
            if d["tag"][j]=="far":
                assert int(code)==255
                continue
            assert str(d["tag"][j])=="fleetfail"+str(code)
            idx=int(d["soff"][j]); count=int(d["nstruct"][j])
            later[int(p)].append((int(i),int(code),d["lo"][j],d["hi"][j],d["shifts"][idx:idx+count],d["lam"][idx:idx+count]))
original={}
for f in results.glob("*/part*.json"):
    doc=json.loads(f.read_bytes())
    for x in doc["fails"]:
        key=doc["part"],x["pos"]
        assert key not in original
        original[key]=x
seen=set()
counts=defaultdict(int)
for p,items in sorted(later.items()):
    codes=np.load(cert/f"part{p:04}.npy")
    meta=json.loads((cert/f"part{p:04}.json").read_bytes())
    prefix=np.r_[1,1+np.cumsum(np.where(codes<11,1,-1))]
    with np.load(cert/f"part{p:04}_multi.npz") as z:
        positions=z["pos"];offset=np.r_[0,np.cumsum(z["nstruct"])]
        shifts,weights=z["shifts"],z["lam"]
        for pos,code,lo,hi,sh,lam in items:
            key=p,pos
            assert key not in seen and key in original
            seen.add(key)
            assert code==codes[pos]==original[key]["code"]
            c,h=leaf_box(prefix,codes,p,pos,meta)
            left=[Fraction(i)-Fraction(j) for i,j in zip(c,h)]
            right=[Fraction(i)+Fraction(j) for i,j in zip(c,h)]
            assert left==list(map(Fraction,lo))==list(map(Fraction,original[key]["lo"]))
            assert right==list(map(Fraction,hi))==list(map(Fraction,original[key]["hi"]))
            if code==253:
                rank=int(np.searchsorted(positions,pos))
                assert positions[rank]==pos
                sl=slice(int(offset[rank]),int(offset[rank+1]))
                assert np.array_equal(shifts[sl],sh) and np.array_equal(weights[sl],lam)
            else:
                assert len(sh)==len(lam)==0
            counts[code]+=1
assert seen==set(original) and len(seen)==2206
print(json.dumps({"pass1_unresolved":len(original),"later_distinct_boxes":len(seen),
    "parts":len(later),"by_code":counts,"all_exact_endpoints_match":True,
    "all_combinations_match_original":True},indent=2))
