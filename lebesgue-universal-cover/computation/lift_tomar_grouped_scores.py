"""Retain high recorded area floors, then check every other closed D3 box at a higher rational target."""
import argparse
import ctypes
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from check_tomar_grouped_residues import sha, scaled, TABLE, FINE, SCALE
from tomar_d3_mask import leaf_mask
from verify_tomar_fan_table import PINNED_MANIFEST, verified_table

p = argparse.ArgumentParser(description=__doc__)
for key in ("ledger", "native", "certdir", "table", "output", "cache"):
    p.add_argument(key, type=Path)
p.add_argument("--target", default="0.8372")
p.add_argument("--parts")
p.add_argument("--calls", type=int, default=4095)
p.add_argument("--depth", type=int, default=30)
a = p.parse_args()
target = F(a.target)
assert F("0.8371") < target < 1
goal = 2*SCALE*SCALE*target
integer_goal = -(-goal.numerator//goal.denominator)
required = target * 10**7
threshold = -(-required.numerator//required.denominator)-8371000+1
ledger = json.loads(a.ledger.read_bytes())
assert ledger["target"] == "0.8371" and ledger["historical_area_inputs"] == [] and sha(a.table) == TABLE
manifest = a.certdir/"SHA256SUMS"
assert sha(manifest) == PINNED_MANIFEST
pinned = {n.removeprefix("./"): h for h,n in (row.split() for row in manifest.read_text().splitlines())}
first = json.loads((a.certdir/"part0000.json").read_bytes())
data, counts = verified_table(a.table, a.certdir, first)
base = scaled(first["hi"])
selected = {r[0]:r[1:] for r in ledger["receipt_hashes"]}
keys = sorted(selected) if a.parts is None else list(map(int, a.parts.split(",")))
assert set(keys) <= set(selected) and len(set(keys)) == len(keys)
a.output.mkdir(parents=True, exist_ok=True)
source = HERE/"tomar_grouped_fan.cpp"
with tempfile.TemporaryDirectory(dir=a.output) as temp:
    binary = Path(temp)/"libtomar_grouped_lift.dylib"
    subprocess.run(["/usr/bin/clang++","-O3","-ftrapv","-std=c++17",f"-fmodules-cache-path={a.cache}",
                    "-shared","-fPIC","-Wl,-install_name,@rpath/libtomar_grouped_lift.dylib",source,"-o",binary],check=True)
    lib = ctypes.CDLL(str(binary))
    vec = lambda t:np.ctypeslib.ndpointer(dtype=t,ndim=1,flags="C_CONTIGUOUS")
    lib.adaptive_trim.argtypes = [vec(np.uint8),ctypes.c_int,vec(np.int64),vec(np.int64),vec(np.int64),
        vec(np.uint8),vec(np.uint8),vec(np.int64),vec(np.int32),vec(np.int32),ctypes.c_int64,ctypes.c_int,ctypes.c_int]
    commitment = {"target":a.target,"minimum_reused_score":threshold,"input_audit_sha256":sha(a.ledger),
                  "code_sha256":sha(source),"driver_sha256":sha(Path(__file__)),"library_sha256":sha(binary),
                  "slice_selector_sha256":sha(HERE/"tomar_d3_mask.py"),"table_sha256":TABLE,
                  "max_depth":a.depth,"max_calls":a.calls}
    for part in keys:
        stem=f"part{part:04}"
        parent=a.native/(stem+".json")
        npz=parent.with_suffix(".npz")
        assert [sha(parent),sha(npz)] == selected[part][:2]
        previous=json.loads(parent.read_bytes())
        header=a.certdir/parent.name
        tree=a.certdir/(stem+".npy")
        assert sha(header)==pinned[header.name]==previous["header_sha256"] and sha(tree)==pinned[tree.name]==previous["codes_sha256"]
        meta=json.loads(header.read_bytes())
        assert meta["part"]==part and meta["hi"]==first["hi"] and meta["lo"]==first["lo"]
        codes=np.load(tree,allow_pickle=False)
        mask=leaf_mask(codes,part,F(meta["hi"][2]))
        assert hashlib.sha256(mask.tobytes()).hexdigest()==previous["slice_mask_sha256"]
        ids=np.flatnonzero((mask==1)|(mask==2))
        with np.load(npz,allow_pickle=False) as saved:
            scores=saved["scores"]
        assert scores.dtype==np.uint16 and len(ids)==len(scores)==previous["all_original"]
        needs=np.ones(len(ids),dtype=bool) if threshold>65535 else scores<threshold
        old=np.zeros_like(mask)
        old[ids[needs]]=1
        lo,hi=[],[]
        for x,y,bit in zip(meta["lo"],meta["hi"],np.unravel_index(part,(2,)*11)):
            l,h=F(x)*FINE,F(y)*FINE
            l,h=((l+h)/2,h) if bit else (l,(l+h)/2)
            assert l.denominator==h.denominator==1
            lo.append(int(l));hi.append(int(h))
        out=np.zeros_like(mask)
        stats=np.zeros(6,dtype=np.int64)
        err=lib.adaptive_trim(codes,len(codes),np.array(lo),np.array(hi),base,old,out,stats,
                              data.ravel(),counts,integer_goal,a.depth,a.calls)
        assert err==0 and np.all(out<=old)
        remaining=np.flatnonzero(out).astype(np.int32)
        assert int(stats[2]+stats[3])==int(old.sum()) and int(stats[3])==len(remaining)
        dest=a.output/parent.name
        objects=dest.with_suffix(".npz")
        assert not dest.exists() and not objects.exists()
        tmp=dest.with_suffix(".tmp.npz")
        np.savez_compressed(tmp,pos=remaining)
        tmp.replace(objects)
        result={**commitment,"part":part,"parent_sha256":selected[part][0],"scores_sha256":selected[part][1],
                "canonical":len(ids),"reused":int((~needs).sum()),"fresh_closed":int(stats[2]),
                "remaining":len(remaining),"checked_boxes":int(stats[0]),"max_depth_used":int(stats[1]),
                "remaining_sha256":sha(objects)}
        tmp=dest.with_suffix(".partial")
        tmp.write_text(json.dumps(result,separators=(",",":"))+"\n")
        tmp.replace(dest)
        print(part,result["canonical"],result["reused"],result["fresh_closed"],result["remaining"],flush=True)
