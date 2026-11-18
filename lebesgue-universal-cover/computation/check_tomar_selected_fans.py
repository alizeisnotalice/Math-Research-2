"""Freshly check selected original boxes against a specified rational target and kernel."""
import argparse, ctypes, hashlib, json, subprocess, sys, time
from fractions import Fraction
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
from tomar_leaf_box import leaf_box
from verify_tomar_fan_table import verified_table, PINNED_MANIFEST
p=argparse.ArgumentParser(description=__doc__)
for k in ["certdir","table","kernel","output","target"]:p.add_argument(k)
p.add_argument("--selected",type=Path,help="conditional historical selection with a weak list")
p.add_argument("--stage",type=Path,help="original parent native-stage JSON/NPZ directory")
p.add_argument("--parts",default=None)
p.add_argument("--calls",type=int,default=32767);p.add_argument("--depth",type=int,default=45)
p.add_argument("--limit",type=int)
a=p.parse_args();cert=Path(a.certdir);out=Path(a.output);out.mkdir(exist_ok=True)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert bool(a.stage)!=bool(a.selected)
source=Path(a.kernel).resolve();assert source.name in ("tomar_curved_fan.cpp","tomar_parabolic_fan.cpp")
library=out/("lib"+source.stem+".dylib")
subprocess.run(["/usr/bin/clang++","-O3","-ftrapv","-std=c++17","-fmodules-cache-path=/Volumes/data/CodexTemp/cache/clang","-shared","-fPIC","-Wl,-install_name,@rpath/"+library.name,source,"-o",library],check=True)
lib=ctypes.CDLL(str(library.resolve()));vec=lambda t:np.ctypeslib.ndpointer(dtype=t,ndim=1,flags="C_CONTIGUOUS")
lib.fan_probe.argtypes=[ctypes.c_int,vec(np.int64),vec(np.int64),vec(np.int64),vec(np.int32),vec(np.int32),ctypes.c_int64,ctypes.c_int,ctypes.c_int,vec(np.int64)]
first=json.loads((cert/"part0000.json").read_bytes());table,counts=verified_table(Path(a.table),cert,first)
manifest=cert/"SHA256SUMS";assert sha(manifest)==PINNED_MANIFEST
expected={f.removeprefix("./"):h for h,f in (line.split() for line in manifest.read_text().splitlines())}
from collections import defaultdict
def scaled(xs):
 v=[Fraction(x)*2**44 for x in xs];assert all(x.denominator==1 for x in v)
 return np.array(list(map(int,v)),dtype=np.int64)
hi0=scaled(first["hi"]);target=2*Fraction(a.target)*(2**28)**2;goal=-(-target.numerator//target.denominator)
rows=defaultdict(list);sources={}
if a.selected:
 selection=json.loads(a.selected.read_bytes());assert selection["target"]==a.target
 for part,pos,code,cls,*_ in selection["weak"]:rows[part].append(pos)
 sources={i:sha(a.selected) for i in rows}
else:
 for f in a.stage.glob("part*.json"):
  rec=json.loads(f.read_bytes());part=rec["part"];assert rec["target"]==a.target and rec["tested"]==rec["all_original"]
  raw=f.with_suffix(".npz");assert sha(raw)==rec["remaining_sha256"]
  with np.load(raw,allow_pickle=False) as d:pos=d["pos"]
  assert len(pos)==rec["remaining"] and np.array_equal(pos,np.unique(pos))
  rows[part]=pos.tolist();sources[part]=sha(f)
chosen=sorted(rows) if a.parts is None else list(map(int,a.parts.split(",")))
for part in chosen:
 start=time.monotonic();stem=f"part{part:04}";header=cert/(stem+".json");codefile=cert/(stem+".npy")
 assert sha(header)==expected[header.name] and sha(codefile)==expected[codefile.name]
 meta=json.loads(header.read_bytes());assert meta["part"]==part and meta["lo"]==first["lo"] and meta["hi"]==first["hi"]
 codes=np.load(codefile,allow_pickle=False);prefix=np.r_[1,1+np.cumsum(np.where(codes<11,1,-1))]
 pos=rows[part];assert len(pos)==len(set(pos))
 if a.limit is not None:pos=pos[:a.limit]
 answers=[]
 for off in range(0,len(pos),32):
  batch=pos[off:off+32];lo=[];hi=[]
  for i in batch:
   assert 0<=i<len(codes) and codes[i]>=253
   c,h=leaf_box(prefix,codes,part,i,meta)
   lo.append(scaled([Fraction(x)-Fraction(y) for x,y in zip(c,h)]));hi.append(scaled([Fraction(x)+Fraction(y) for x,y in zip(c,h)]))
  stats=np.zeros((len(batch),4),dtype=np.int64)
  assert lib.fan_probe(len(batch),np.array(lo).ravel(),np.array(hi).ravel(),hi0,table.ravel(),counts,goal,a.depth,a.calls,stats.ravel())==0
  answers.extend([[i,int(codes[i]),*map(int,res)] for i,res in zip(batch,stats)])
 record={"part":part,"target":a.target,"scope":"selected raw leaf indices","selection_sha256":sources[part],"source_sha256":sha(source),
  "driver_sha256":sha(__file__),"decoder_sha256":sha(Path(__file__).with_name("tomar_leaf_box.py")),"table_sha256":sha(a.table),
  "library_sha256":sha(library),"header_sha256":sha(header),"codes_sha256":sha(codefile),"depth":a.depth,"calls":a.calls,
  "checked":len(pos),"closed":sum(row[-1] for row in answers),"seconds":round(time.monotonic()-start,3),"leaves":answers}
 dest=out/(stem+".json");tmp=dest.with_suffix(".partial");tmp.write_text(json.dumps(record,separators=(",",":"))+"\n");tmp.replace(dest)
 print(part,len(pos),record["closed"],record["seconds"],flush=True)
