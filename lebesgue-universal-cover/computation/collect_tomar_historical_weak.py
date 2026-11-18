"""Identify original or D3-slice weak leaves relative to Tomar's frozen historical arithmetic."""
import argparse
from collections import Counter
from fractions import Fraction
import hashlib,json,subprocess,tarfile
from pathlib import Path
import numpy as np
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
from tomar_five_orbit_slice import leaf_mask, CUT
from flint import arb,ctx
p=argparse.ArgumentParser(description=__doc__)
for name in ("checkout","certdir","archive","results","masks","mask_index","out"):p.add_argument(name,type=Path)
p.add_argument("--target",action="append",required=True)
p.add_argument("--recompute-slice",action="store_true")
a=p.parse_args();a.out.mkdir(exist_ok=True)
def digest(path):
 with path.open("rb") as stream:return hashlib.file_digest(stream,"sha256").hexdigest()
revision="10af98606c98ca9ffd0d9f402cae1e2efe74b792"
assert digest(a.archive)=="0cf9696eb47681382492bc547d44519beeb69f492ee0cf26c520ec8254e1562b"
manifest=a.certdir/"SHA256SUMS"
assert digest(manifest)=="56cbe5dfdadc6b5aa07cf6ec9acf4f637a3c4fa78db54df0ae2d78916d292410"
pinned={file.removeprefix("./"):h for h,file in (x.split() for x in manifest.read_text().splitlines())}
files={}
with tarfile.open(a.archive,"r|xz") as archive:
 for member in archive:
  if member.isfile() and member.name.endswith((".json","_margins.npy")) and "/part" in member.name:
   name=member.name.removeprefix("fleet11/results/")
   path=a.results/name
   assert digest(path)==hashlib.file_digest(archive.extractfile(member),"sha256").hexdigest()
   files[name]=path
assert digest(a.mask_index)=="ec483a5d64c86fc1d305c63081e08d66fd128bb2033622764ee26b2d45377e02"
ctx.prec=128
assert arb(CUT.numerator)/CUT.denominator>arb.pi()/15
index=json.loads(a.mask_index.read_bytes());assert index["manifest"]==digest(manifest) and index["active_roots"]==1024
assert index["counts"]==[4002817,2707497,9483088] and index["cut"]=="1308996939/6250000000"
later={}
root=a.checkout/"independent/full11/pass23"
def frozen(path):
 assert subprocess.check_output(["git","-C",str(a.checkout),"show",revision+":"+str(path.relative_to(a.checkout))])==path.read_bytes()
for file in sorted((root/"results").glob("*.jsonl")):
 frozen(file)
 side=root/("fails_"+file.stem.split("_",1)[1]+".npz");frozen(side)
 with np.load(side,allow_pickle=False) as z:
  for line in file.read_bytes().splitlines():
   rec=json.loads(line);ix=rec["idx"]
   key=(int(z["part"][ix]),int(z["pos"][ix]),int(z["code"][ix]))
   assert key[:2]==(rec["part"],rec["pos"])
   if rec["ok"]:later[key]=max(later.get(key,0),float(rec["margin"]))
targets={s:Fraction(s)-Fraction(8366,10000) for s in a.target};assert all(gap>0 for gap in targets.values())
weak={s:[] for s in targets};counts={s:Counter() for s in targets};total=0;scope=0;seen=set();used=set()
for fn,f in sorted(files.items()):
 if not fn.endswith(".json"):continue
 rec=json.loads(f.read_bytes());part=rec["part"];assert part not in seen;seen.add(part)
 val=files[fn[:-5]+"_margins.npy"];v=np.load(val,allow_pickle=False)
 codesfile=a.certdir/f"part{part:04}.npy";assert digest(codesfile)==pinned[codesfile.name]
 code=np.load(codesfile,allow_pickle=False);pos=np.flatnonzero(code>=253)
 assert len(v)==len(pos) and v.dtype==np.float64
 mfile=a.masks/f"part{part:04}.npy";row=index["rows"][part];mask=np.load(mfile,allow_pickle=False);assert mask.dtype==np.uint8 and mask.shape==code.shape
 assert row[0]==part and hashlib.sha256(mask.tobytes()).hexdigest()==row[-1]
 if a.recompute_slice:
  head=a.certdir/f"part{part:04}.json";assert digest(head)==pinned[head.name]
  assert np.array_equal(mask,leaf_mask(code,part,Fraction(json.loads(head.read_bytes())["hi"][2])))
 m=mask[pos];assert [int((m==k).sum()) for k in (1,2,3)]==row[1:4]
 fail=np.flatnonzero(np.isnan(v));later_rank=set(map(int,fail))
 for rank in fail:
  i=int(pos[rank]);key=(part,i,int(code[i]));assert key in later;v[rank]=later[key];used.add(key)
 assert np.isfinite(v).all() and (v>0).all()
 low=np.nextafter(v,-np.inf);total+=len(v);scope+=int(np.count_nonzero(m<3))
 for target,gap in targets.items():
  for j in np.flatnonzero(low<=np.nextafter(float(gap),np.inf)):
   if Fraction(float(low[j]))>gap:continue
   cls=int(m[j]);source="later" if int(j) in later_rank else "pass1"
   counts[target][source+"_weak_full"]+=1;counts[target][f"category_{cls}"]+=1
   if cls<=2:
    i=int(pos[j]);weak[target].append([part,i,int(code[i]),cls,source,float(v[j]).hex()])
assert seen==set(range(2048)) and total==16193402 and scope==6710314 and len(used)==2206
for t in targets:
 record={"status":"strong cases conditional on Tomar's pass1/pass23 area arithmetic","target":t,
  "historic_tau":"0.8366","full_leaves":total,"canonical_outer_leaves":scope,"data":counts[t],
  "canonical_weak":len(weak[t]),"canonical_strong":scope-len(weak[t]),"pass1_archive_sha256":digest(a.archive),
  "input_manifest_sha256":digest(manifest),"slice_index_sha256":digest(a.mask_index),"upstream_revision":revision,
  "columns":["part","pos","code","slice","history_pass","historical_margin_hex"],"weak":weak[t]}
 dest=a.out/("weak-"+t.replace(".","")+".json");dest.write_text(json.dumps(record,separators=(",",":"))+"\n")
 print(t,dict(counts[t]),len(weak[t]),digest(dest),flush=True)
