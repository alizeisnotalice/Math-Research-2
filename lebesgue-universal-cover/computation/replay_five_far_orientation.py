"""Audit or independently recompute all 155 original leaves in the far ninth-angle component."""
import argparse
from fractions import Fraction
import os
import platform
import shutil
from collections import defaultdict
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from zipfile import ZipFile

for name in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","VECLIB_MAXIMUM_THREADS"):
    os.environ[name]="1"
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from tomar_leaf_box import leaf_box
import locate_tomar_five_basin as localize
def digest(b):
    return sha256(b).hexdigest()

def main():
    if not __debug__:raise SystemExit("Assertions required.")
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("checkout","cert","high","archive","scratch","cache"):p.add_argument(key,type=Path)
    p.add_argument("--replay",action="store_true")
    a=p.parse_args()
    assert digest(a.archive.read_bytes())=="1600ee09d8d929e4286b0e116c66528efc0c413d5acfef7fc504dfe5a29d5342"
    index=HERE/"five-orientation-branches.json"
    raw=subprocess.check_output([sys.executable,"-I","-B",HERE/"certify_five_orientation_branches.py",a.cert,a.high])
    assert raw==index.read_bytes() and digest(raw)=="654c72a37c28239680df4880a96319e23ed817ff9cc039b58f644bfc27935f2f"
    ix=json.loads(raw);target="0.83729500622102525"
    far={tuple(k) for group in ix["components"][1:] for k in group};assert len(far)==155
    parts=defaultdict(list)
    for k,i in sorted(far):parts[k].append(i)
    runner=HERE/"locate_tomar_five_basin.py";native=HERE/"tomar_dynamic_affine.cpp";driver=HERE/"tomar_dynamic_five.py"
    sum_vis=0;closed=defaultdict(int);replayed=0
    with ZipFile(a.archive)as z,tempfile.TemporaryDirectory(dir=a.scratch)as tmp:
        metadata=json.loads((a.high/"five-body-contact-basin.json").read_bytes())
        basin=localize.inward_box(metadata)
        if a.replay:
            rev="10af98606c98ca9ffd0d9f402cae1e2efe74b792";up=a.checkout.resolve()
            assert subprocess.check_output(["git","-C",up,"rev-parse","HEAD"],text=True).strip()==rev
            for name in ("gen8.py","gen11.py"):
                assert (up/"independent"/name).read_bytes()==subprocess.check_output(["git","-C",up,"show",rev+":independent/"+name])
            lib=Path(tmp)/"libtomar_unique.dylib"
            compiler=shutil.which("clang++")
            assert compiler
            flags=["-shared","-fPIC"]
            if platform.system()=="Darwin":flags+=["-Wl,-install_name,@rpath/libtomar_unique.dylib"]
            subprocess.run([compiler,"-O3","-ftrapv","-std=c++17",f"-fmodules-cache-path={a.cache}",*flags,native,"-o",lib],check=True)
            sys.path.insert(0,str(up/"independent"));import gen11
            localize.dynamic.configure(up,a.cert,lib,target);assert localize.dynamic.gen8.placements is gen11.placements
        assert sorted(z.namelist())==sorted(f"part{k:04}-pos{i}.json"for k,i in far)
        for k,ids in parts.items():
            header=a.cert/f"part{k:04}.json";codes=np.load(header.with_suffix(".npy"),allow_pickle=False)
            meta=json.loads(header.read_bytes());prefix=np.r_[1,1+np.cumsum(np.where(codes<11,1,-1))]
            for pos in ids:
                name=f"part{k:04}-pos{pos}.json";b=z.read(name);x=json.loads(b)
                assert (x["part"],x["pos"])==(k,pos) and x["target"]==target and x["budget"]==4096
                assert x["parent_sha256"]==ix["parent_sha256"] and x["tree_manifest_sha256"]==ix["input_manifest_sha256"]
                assert x["source_sha256"]==digest(native.read_bytes()) and x["dynamic_sha256"]==digest(driver.read_bytes())
                assert x["driver_sha256"]==digest(runner.read_bytes()) and x["basin_sha256"]==ix["basin_sha256"]
                c,w=leaf_box(prefix,codes,k,pos,meta)
                root=[(Fraction(v)+s*Fraction(d))*2**44 for s in (-1,1)for v,d in zip(c,w)]
                assert all(n.denominator==1 for n in root) and x["root"]==root
                q=x["result"];assert q["live"]==[] and not set(q["closed"])-{"area1","area2","area3","radial"}
                ends=sum(q["closed"].values());assert ends==q["splits"]+1
                assert q["visits"]==ends+q["splits"]+q["contractions"]
                sum_vis+=q["visits"]
                for label,n in q["closed"].items():closed[label]+=n
                if a.replay:
                    actual=localize.cover(x["root"],basin,4096)
                    assert actual==q;replayed+=1
    print(json.dumps(dict(target=target,original_leaves=len(far),canonical_far_angle_strict_upper="9/200",
        excluded_all_far=True,visits=sum_vis,terminal_counts=closed,trees_replayed_in_this_run=replayed,
        archive_sha256=digest(a.archive.read_bytes()),archive_bytes=a.archive.stat().st_size),indent=2))
if __name__=="__main__":main()
