"""Verify the pointwise stronger-score partition against an authenticated grouped-fan parent ledger."""
import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
from flint import arb,ctx

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from check_tomar_grouped_residues import sha, TABLE
from tomar_d3_mask import leaf_mask

p = argparse.ArgumentParser(description=__doc__)
for key in ("parent_ledger","parent_native","certdir","higher","scratch","cache"):
    p.add_argument(key,type=Path)
p.add_argument("--target",required=True)
p.add_argument("--calls",type=int,default=1048575)
p.add_argument("--depth",type=int,default=60)
p.add_argument("--partial",action="store_true")
a = p.parse_args()
source = HERE/"tomar_grouped_fan.cpp"
with tempfile.TemporaryDirectory(dir=a.scratch) as tmp:
    library=Path(tmp)/"libtomar_grouped_lift.dylib"
    subprocess.run(["/usr/bin/clang++","-O3","-ftrapv","-std=c++17",f"-fmodules-cache-path={a.cache}",
        "-shared","-fPIC","-Wl,-install_name,@rpath/libtomar_grouped_lift.dylib",source,"-o",library],check=True)
    binary_hash=sha(library)
parents=json.loads(a.parent_ledger.read_bytes())
assert parents["target"]=="0.8371" and parents["historical_area_inputs"]==[]
target=Fraction(a.target);scaled=target*10**7
assert target>Fraction("0.8371")
ctx.prec=128
assert len(parents["root_exterior"])==4 and all(arb(x)>arb(a.target) for x in parents["root_exterior"])
threshold=-(-scaled.numerator//scaled.denominator)-8371000+1
rows={r[0]:r[1:] for r in parents["receipt_hashes"]}
parts=sorted(int(f.stem[4:]) for f in a.higher.glob("part*.json"))
assert parts and set(parts)<=set(rows)
if not a.partial:
    assert parents["complete"] and len(rows)==len(parts)==1024
totals=Counter();committed=[]
for part in parts:
    stem=f"part{part:04}"
    prior=a.parent_native/(stem+".json");values=prior.with_suffix(".npz")
    assert sha(prior)==rows[part][0] and sha(values)==rows[part][1]
    low=json.loads(prior.read_bytes())
    tree=a.certdir/(stem+".npy");header=a.certdir/prior.name
    assert sha(tree)==low["codes_sha256"] and sha(header)==low["header_sha256"]
    meta=json.loads(header.read_bytes());codes=np.load(tree,allow_pickle=False)
    cls=leaf_mask(codes,part,Fraction(meta["hi"][2]))
    assert hashlib.sha256(cls.tobytes()).hexdigest()==low["slice_mask_sha256"]
    pos=np.flatnonzero((cls==1)|(cls==2))
    with np.load(values,allow_pickle=False) as old:
        score=old["scores"]
    assert score.dtype==np.uint16 and score.shape==pos.shape and len(pos)==low["all_original"]
    strong=np.zeros(len(pos),dtype=bool) if threshold>65535 else score>=threshold
    f=a.higher/prior.name;v=f.with_suffix(".npz");r=json.loads(f.read_bytes())
    for key,expected in {"target":a.target,"minimum_reused_score":threshold,"input_audit_sha256":sha(a.parent_ledger),
            "code_sha256":sha(source),"driver_sha256":sha(HERE/"lift_tomar_grouped_scores.py"),
            "library_sha256":binary_hash,"slice_selector_sha256":sha(HERE/"tomar_d3_mask.py"),"table_sha256":TABLE,
            "max_depth":a.depth,"max_calls":a.calls,"part":part,"parent_sha256":rows[part][0],
            "scores_sha256":rows[part][1],"canonical":len(pos),"reused":int(strong.sum()),"remaining_sha256":sha(v)}.items():
        assert r[key]==expected,(part,key)
    with np.load(v,allow_pickle=False) as z:
        assert z.files==["pos"];left=z["pos"]
    assert left.dtype==np.int32 and np.array_equal(left,np.unique(left)) and np.isin(left,pos[~strong]).all()
    assert len(left)==r["remaining"] and len(pos)==r["reused"]+r["fresh_closed"]+len(left)
    assert 0<=r["max_depth_used"]<=a.depth and 0<=r["checked_boxes"]
    totals.update({k:r[k] for k in ("canonical","reused","fresh_closed","remaining","checked_boxes")})
    committed.append([part,sha(f),sha(v)])
if not a.partial:
    assert totals["canonical"]==6710314 and totals["remaining"]==0
print(json.dumps({"target":a.target,"authenticated_parent_ledger_sha256":sha(a.parent_ledger),
   "historical_area_inputs":[],"complete":not a.partial,"parts":len(parts),"totals":dict(totals),
   "higher_receipt_hashes":committed},separators=(",",":")))
