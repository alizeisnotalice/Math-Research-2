"""Recompute whole original five-body boxes by dynamic affine lower forms at the 0.83725 target."""
import argparse
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import tempfile

for key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","VECLIB_MAXIMUM_THREADS"):
    os.environ[key]="1"
import numpy as np
import scipy
import flint

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from tomar_d3_mask import leaf_mask
from tomar_leaf_box import leaf_box
from verify_tomar_fan_table import PINNED_MANIFEST
import tomar_dynamic_five as dynamic

def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("checkout","cert","slice_index","selection","output","cache"):p.add_argument(key,type=Path)
    p.add_argument("--parts")
    a=p.parse_args()
    sel=json.loads(a.selection.read_bytes());target="0.83725"
    assert Fraction(sel["target"])>Fraction(target) and sel["totals"][-2]==sum(map(len,sel["retained"].values()))
    assert sha(a.slice_index)=="ec483a5d64c86fc1d305c63081e08d66fd128bb2033622764ee26b2d45377e02"
    index=json.loads(a.slice_index.read_bytes())
    manifest=a.cert/"SHA256SUMS";assert sha(manifest)==PINNED_MANIFEST
    pinned={n.removeprefix("./"):h for h,n in (s.split()for s in manifest.read_text().splitlines())}
    checkout=a.checkout.resolve()
    rev=subprocess.check_output(["git","-C",checkout,"rev-parse","HEAD"],text=True).strip()
    assert rev=="10af98606c98ca9ffd0d9f402cae1e2efe74b792"
    proposal=checkout/"independent/gen8.py"
    official=subprocess.check_output(["git","-C",checkout,"show",rev+":independent/gen8.py"])
    assert hashlib.sha256(official).hexdigest()==sha(proposal)
    a.output.mkdir(parents=True,exist_ok=True)
    source=HERE/"tomar_dynamic_affine.cpp"
    with tempfile.TemporaryDirectory(dir=a.output) as tmp:
        lib=Path(tmp)/"libtomar_dynamic_five.dylib"
        subprocess.run(["/usr/bin/clang++","-O3","-ftrapv","-std=c++17",f"-fmodules-cache-path={a.cache}",
            "-shared","-fPIC","-Wl,-install_name,@rpath/libtomar_dynamic_five.dylib",source,"-o",lib],check=True)
        dynamic.configure(checkout,a.cert,lib,target)
        shared={"target":target,"selection_sha256":sha(a.selection),"slice_index_sha256":sha(a.slice_index),
            "input_manifest_sha256":PINNED_MANIFEST,"source_sha256":sha(source),"library_sha256":sha(lib),
            "affine_driver_sha256":sha(HERE/"tomar_dynamic_five.py"),"driver_sha256":sha(Path(__file__)),
            "decoder_sha256":sha(HERE/"tomar_leaf_box.py"),"proposal_sha256":sha(proposal),"upstream":rev,
            "flint":flint.__version__,"numpy":np.__version__,"scipy":scipy.__version__,"python":platform.python_version(),
            "moves_seed":20261118,"max_calls":8192,"max_depth":90,"moves":65}
        parts=sorted(map(int,sel["retained"])) if not a.parts else sorted(map(int,a.parts.split(",")))
        assert all(str(k)in sel["retained"]for k in parts)
        for part in parts:
            stem=f"part{part:04}";h=a.cert/(stem+".json");tree=h.with_suffix(".npy")
            assert sha(h)==pinned[h.name] and sha(tree)==pinned[tree.name]
            head=json.loads(h.read_bytes());codes=np.load(tree,allow_pickle=False)
            labels=leaf_mask(codes,part,Fraction(head["hi"][2]))
            assert hashlib.sha256(labels.tobytes()).hexdigest()==index["rows"][part][-1]
            prefix=np.r_[1,1+np.cumsum(np.where(codes<11,1,-1))]
            pos=sel["retained"][str(part)]
            assert pos==sorted(set(pos)) and np.isin(labels[pos],(1,2)).all()
            rows=[]
            for ident in pos:
                c,w=leaf_box(prefix,codes,part,ident,head)
                v=[(Fraction(x)-Fraction(y))*(2**44)for x,y in zip(c,w)]+[(Fraction(x)+Fraction(y))*(2**44)for x,y in zip(c,w)]
                assert all(x.denominator==1 for x in v)
                q=dynamic.work((part,ident,list(map(int,v)),8192,90))
                rows.append(q)
                print(part,ident,q["vis"],len(q["remaining"]),flush=True)
            receipt={**shared,"part":part,"header_sha256":sha(h),"codes_sha256":sha(tree),"rows":rows}
            out=a.output/(stem+".json");tmpout=out.with_suffix(".partial")
            tmpout.write_text(json.dumps(receipt,separators=(",",":"))+"\n");tmpout.replace(out)

if __name__=="__main__":main()
