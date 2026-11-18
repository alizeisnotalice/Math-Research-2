"""Search the strict five-body minimizer on the 501 original boxes using certified exclusions."""
import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile

for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[name] = "1"
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from tomar_leaf_box import leaf_box
from verify_tomar_fan_table import PINNED_MANIFEST
import tomar_dynamic_five as dynamic

S = 1 << 44
TARGET = "0.83729500622102525"
WEIGHTS = np.array([1, 1, .16, 1, 1, .11, 1, 1, .09, 1, 1])
RADII = [1202425916136, 1169880371954, 1056234850107, 1002578682672]
def digest(file):
    return sha256(file.read_bytes()).hexdigest()

def inward_box(basin):
    cent = list(map(Fraction, basin["raw_center"]))
    width = list(map(Fraction, basin["raw_halfwidth"]))
    return (np.array([math.ceil((c-h)*S) for c,h in zip(cent,width)],dtype=np.int64),
            np.array([math.floor((c+h)*S) for c,h in zip(cent,width)],dtype=np.int64))

def clip_domain(l,h):
    h[2] = min(h[2],3684498829194)
    if l[2]>h[2]:return False
    for j,r in enumerate(RADII):
        i=3*j
        x=max(0,int(l[i]),-int(h[i]));y=max(0,int(l[i+1]),-int(h[i+1]))
        if x*x+y*y>r*r:return False
        for k,d in ((i,y),(i+1,x)):
            z=r*r-d*d
            m=math.isqrt(z);cap=m+int(m*m<z)
            l[k]=max(l[k],-cap);h[k]=min(h[k],cap)
    return True

def cover(v,basin,budget):
    low,high=basin
    l=np.array(v[:11],dtype=np.int64);h=np.array(v[11:],dtype=np.int64)
    pending=[(l,h,0)]
    visits=splits=shortened=0;counts=Counter();fixed=[]
    while pending and visits<budget:
        l,h,depth=pending.pop();visits+=1
        if not clip_domain(l,h):
            counts["radial"]+=1;continue
        if np.all(low<=l) and np.all(h<=high):
            counts["basin"]+=1;continue
        ok,axis,L,H,_=dynamic.contract(l,h,65)
        if ok:
            counts[f"area{ok}"]+=1;continue
        if np.all(low<=L) and np.all(H<=high):
            counts["basin"]+=1;continue
        width=(h-l)*WEIGHTS;eff=(H-L)*WEIGHTS
        if (width-eff).max()>.1*max(1,width.max()):
            shortened+=1;pending.append((L,H,depth+1));continue
        if eff[axis]*24<eff.max() or H[axis]-L[axis]<2:
            possible=np.where(H-L>=2,eff,-1)
            axis=int(np.argmax(possible))
        if H[axis]-L[axis]<2:
            fixed.append((L,H,depth));continue
        mid=L[axis]+(H[axis]-L[axis])//2
        C=L.copy();C[axis]=mid;B=H.copy();B[axis]=mid
        splits+=1;pending.extend([(C,H,depth+1),(L,B,depth+1)])
    live=[[*l.tolist(),*h.tolist(),n] for l,h,n in pending+fixed]
    assert sum(counts.values())+len(live)==splits+1
    return dict(visits=visits,splits=splits,contractions=shortened,closed=counts,live=live)

def main():
    if not __debug__:raise SystemExit("Assertions required")
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("checkout","cert","high","output","cache"):p.add_argument(key,type=Path)
    p.add_argument("--parts",required=True);p.add_argument("--budget",type=int,default=2048)
    p.add_argument("--positions",help="selected original indices for a single part")
    a=p.parse_args()
    parent=a.high/"five-global-localization.json"
    basin=a.high/"five-body-contact-basin.json"
    assert digest(parent)=="d85941b19b47cda83fc86268daa50569db2d80118b7fc722d4bf90adfa457504"
    expected_basin=subprocess.check_output(["git","-C",a.high,"show","d8c1b4891effc892a2e08b9568df90a9ddae38c4:lebesgue-universal-cover/computation/five-body-contact-basin.json"])
    assert basin.read_bytes()==expected_basin
    metadata=json.loads(parent.read_bytes())
    assert Fraction(metadata["target"])>Fraction(TARGET)
    box=inward_box(json.loads(expected_basin))
    checkout=a.checkout.resolve();rev="10af98606c98ca9ffd0d9f402cae1e2efe74b792"
    assert subprocess.check_output(["git","-C",checkout,"rev-parse","HEAD"],text=True).strip()==rev
    for f in ("gen8.py","gen11.py"):
        assert (checkout/"independent"/f).read_bytes()==subprocess.check_output(["git","-C",checkout,"show",rev+":independent/"+f])
    m=a.cert/"SHA256SUMS";assert digest(m)==PINNED_MANIFEST
    pins={name.removeprefix("./"):h for h,name in (r.split() for r in m.read_text().splitlines())}
    assert digest(a.cert/"part0000.json")==pins["part0000.json"]
    a.output.mkdir(exist_ok=True,parents=True)
    import platform
    src=HERE/"tomar_dynamic_affine.cpp"
    with tempfile.TemporaryDirectory(dir=a.output) as t:
        lib=Path(t)/"libtomar_unique.dylib"
        subprocess.run(["/usr/bin/clang++","-O3","-ftrapv","-std=c++17",f"-fmodules-cache-path={a.cache}",
            "-shared","-fPIC","-Wl,-install_name,@rpath/libtomar_unique.dylib",src,"-o",lib],check=True)
        sys.path.insert(0,str(checkout/"independent"));import gen11
        dynamic.configure(checkout,a.cert,lib,TARGET);assert dynamic.gen8.placements is gen11.placements
        common={"target":TARGET,"parent_sha256":digest(parent),"basin_sha256":digest(basin),"source_sha256":digest(src),
            "kernel_sha256":digest(lib),"dynamic_sha256":digest(HERE/"tomar_dynamic_five.py"),"driver_sha256":digest(Path(__file__)),
            "tree_manifest_sha256":PINNED_MANIFEST,"upstream":rev,"budget":a.budget,"python":platform.python_version()}
        parts=list(map(int,a.parts.split(",")))
        assert not a.positions or len(parts)==1
        for part in parts:
            ids=metadata["retained"][str(part)]
            if a.positions:ids=list(map(int,a.positions.split(",")));assert set(ids)<=set(metadata["retained"][str(part)])
            header=a.cert/f"part{part:04}.json";tree=header.with_suffix(".npy")
            assert digest(header)==pins[header.name] and digest(tree)==pins[tree.name]
            meta=json.loads(header.read_bytes());codes=np.load(tree,allow_pickle=False)
            pref=np.r_[1,1+np.cumsum(np.where(codes<11,1,-1))]
            for pos in ids:
                dst=a.output/f"part{part:04}-pos{pos}.json";assert not dst.exists()
                c,w=leaf_box(pref,codes,part,pos,meta)
                roots=[(Fraction(x)+sg*Fraction(y))*S for sg in (-1,1) for x,y in zip(c,w)]
                assert all(z.denominator==1 for z in roots)
                root=list(map(int,roots))
                result=cover(root,box,a.budget)
                row={**common,"part":part,"pos":pos,"root":root,"result":result}
                out=dst.with_suffix(".partial");out.write_text(json.dumps(row,separators=(",",":"))+"\n");out.replace(dst)
                print(part,pos,result["visits"],len(result["live"]),result["closed"].get("basin",0),flush=True)
if __name__=="__main__":main()
