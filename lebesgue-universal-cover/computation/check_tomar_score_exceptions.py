"""Independently recover and verify every outstanding original position of the scored threshold lift."""
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

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from check_tomar_grouped_residues import sha, scaled, TABLE, FINE, SCALE
from tomar_d3_mask import leaf_mask
from verify_tomar_fan_table import PINNED_MANIFEST, verified_table

p=argparse.ArgumentParser(description=__doc__)
for key in ("localization","cert","table","slice_index","kernel","output","cache"):
    p.add_argument(key,type=Path)
p.add_argument("--target",default="0.8372")
p.add_argument("--calls",type=int,default=16777216)
p.add_argument("--depth",type=int,default=100)
p.add_argument("--parts")
a=p.parse_args()
target=F(a.target)
parent=json.loads(a.localization.read_bytes())
assert parent["kind"]=="unclosed original leaves of the scored threshold lift"
assert F(parent["target"])==target==F("0.8372")
assert sum(len(pos) for pos in parent["retained"].values())==parent["count"]
assert a.kernel.name=="tomar_fine_angle_fan.cpp"
assert sha(a.table)==TABLE and sha(a.slice_index)=="ec483a5d64c86fc1d305c63081e08d66fd128bb2033622764ee26b2d45377e02"
index=json.loads(a.slice_index.read_bytes())
manifest=a.cert/"SHA256SUMS";assert sha(manifest)==PINNED_MANIFEST
pinned={n.removeprefix("./"):h for h,n in (row.split() for row in manifest.read_text().splitlines())}
base=json.loads((a.cert/"part0000.json").read_bytes())
data,counts=verified_table(a.table,a.cert,base)
base_lo,base_hi=scaled(base["lo"]),scaled(base["hi"])
mids=(base_lo+base_hi)//2
assert np.all(2*mids==base_lo+base_hi)
a.output.mkdir(parents=True,exist_ok=True)
with tempfile.TemporaryDirectory(dir=a.output) as tmp:
    library=Path(tmp)/"libtomar_score_exceptions.dylib"
    subprocess.run(["/usr/bin/clang++","-O3","-ftrapv","-std=c++17",f"-fmodules-cache-path={a.cache}",
        "-shared","-fPIC","-Wl,-install_name,@rpath/libtomar_score_exceptions.dylib",a.kernel,"-o",library],check=True)
    lib=ctypes.CDLL(str(library));vec=lambda t:np.ctypeslib.ndpointer(dtype=t,ndim=1,flags="C_CONTIGUOUS")
    fn=lib.adaptive_trim
    fn.argtypes=[vec(np.uint8),ctypes.c_int,vec(np.int64),vec(np.int64),vec(np.int64),vec(np.uint8),vec(np.uint8),
        vec(np.int64),vec(np.int32),vec(np.int32),ctypes.c_int64,ctypes.c_int,ctypes.c_int]
    threshold=2*target*SCALE**2;goal=-(-threshold.numerator//threshold.denominator)
    keys=sorted(map(int,parent["retained"])) if not a.parts else list(map(int,a.parts.split(",")))
    # Short shards give early feedback on the tail, long ones remain independently repeatable.
    for k in sorted(keys,key=lambda k:len(parent["retained"][str(k)])):
        stem=f"part{k:04}";tree=a.cert/(stem+".npy");header=a.cert/(stem+".json")
        assert sha(tree)==pinned[tree.name] and sha(header)==pinned[header.name]
        meta=json.loads(header.read_bytes());assert meta["part"]==k and meta["hi"]==base["hi"] and meta["lo"]==base["lo"]
        codes=np.load(tree,allow_pickle=False);classes=leaf_mask(codes,k,F(meta["hi"][2]))
        assert hashlib.sha256(classes.tobytes()).hexdigest()==index["rows"][k][-1]
        pos=parent["retained"][str(k)]
        assert pos==sorted(set(pos)) and np.isin(classes[pos],(1,2)).all()
        bits=np.array(np.unravel_index(k,(2,)*11),bool);lo,hi=base_lo.copy(),base_hi.copy()
        lo[bits],hi[~bits]=mids[bits],mids[~bits]
        old=np.zeros_like(classes);old[pos]=1;mask=np.zeros_like(old);st=np.zeros(6,np.int64)
        assert fn(codes,len(codes),lo,hi,base_hi,old,mask,st,data.ravel(),counts,goal,a.depth,a.calls)==0
        left=np.flatnonzero(mask).tolist()
        assert np.all(mask<=old) and len(left)==int(st[3]) and int(st[2])+len(left)==len(pos)
        result={"part":k,"target":a.target,"input_kind":parent["kind"],"selection_sha256":sha(a.localization),
            "slice_index_sha256":sha(a.slice_index),"header_sha256":sha(header),"codes_sha256":sha(tree),
            "source_sha256":sha(a.kernel),"library_sha256":sha(library),"driver_sha256":sha(Path(__file__)),
            "table_sha256":TABLE,"depth":a.depth,"calls":a.calls,"tested":len(pos),"closed":int(st[2]),
            "checked_boxes":int(st[0]),"maximum_depth":int(st[1]),"remaining":left}
        file=a.output/(stem+".json");assert not file.exists()
        temp=file.with_suffix(".partial");temp.write_text(json.dumps(result,separators=(",",":"))+"\n");temp.replace(file)
        print(k,len(pos),result["closed"],result["checked_boxes"],flush=True)
