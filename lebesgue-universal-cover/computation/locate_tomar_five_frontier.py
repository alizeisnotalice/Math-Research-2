"""Recheck a genuine prior frontier using the certified contact basin and exact-arc material chains."""
import argparse
from collections import Counter
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

for key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","VECLIB_MAXIMUM_THREADS"):
    os.environ[key]="1"
import numpy as np
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE))
import locate_tomar_five_basin as prior
from tomar_five_contact_material import Material,CurvedBasin,contract_forms

def digest(p):return sha256(p.read_bytes()).hexdigest()

class Step:
    def __init__(self,source,straight):
        self.mat=Material(source);self.chart=CurvedBasin(source);self.straight=straight
        self.center=np.asarray([round(float(v)*prior.S)for v in self.chart.cb.RAW],np.int64)
    def inside(self,l,h):
        return (np.all(self.straight[0]<=l)and np.all(h<=self.straight[1])) or self.chart.contains(l,h)
    def material(self,l,h):
        # The two guards only choose which inequalities to attempt.
        if np.max(np.maximum(abs(l[:2]-self.center[:2]),abs(h[:2]-self.center[:2])))>prior.S/100:return False,l,h,0
        mats=0
        for label in ("fixed","near"):
            if label=="fixed":a=self.mat.fixed_anchor
            else:
                if not mats and np.max((h-l)*prior.WEIGHTS)>prior.S/1500:break
                try:a=self.mat.near_anchor(l,h)
                except (ValueError,AssertionError,ZeroDivisionError,OverflowError):continue
            f=self.mat.formsets(l,h,a,prior.TARGET)
            if not f:continue
            mats+=1
            ok,l,h=contract_forms(f,l,h)
            if ok:return True,l,h,mats
        return False,l,h,mats
    def tree(self,v,budget):
        stack=[(np.array(v[:11],np.int64),np.array(v[11:22],np.int64),v[22])]
        end=Counter();visits=splits=cuts=trials=basin_faces=0;stalled=[]
        while stack and visits<budget:
            l,h,depth=stack.pop();visits+=1
            if not prior.clip_domain(l,h):end["radial"]+=1;continue
            if self.inside(l,h):end["basin"]+=1;continue
            t,axis,L,H,_=prior.dynamic.contract(l,h,65)
            if t:end[f"area{t}"]+=1;continue
            if self.inside(L,H):end["basin"]+=1;continue
            ok,L,H,m=self.material(L,H);trials+=m
            if ok:end["material"]+=1;continue
            if self.inside(L,H):end["basin"]+=1;continue
            w=(h-l)*prior.WEIGHTS;e=(H-L)*prior.WEIGHTS
            if np.max(w-e)>.1*max(1,np.max(w)):
                cuts+=1;stack.append((L,H,depth+1));continue
            if e[axis]*24<e.max()or H[axis]-L[axis]<2:
                axis=int(np.argmax(np.where(H-L>=2,e,-1)))
            if H[axis]-L[axis]<2:
                stalled.append((L,H,depth));continue
            mid=L[axis]+(H[axis]-L[axis])//2
            bl,bh=self.straight
            if np.all(np.maximum(L,bl)<=np.minimum(H,bh)) and np.max((H-L)*prior.WEIGHTS)<prior.S/1000:
                candidates=[(min(int(v-L[j]),int(H[j]-v))*prior.WEIGHTS[j],j,int(v))
                            for j in range(11) for v in (bl[j],bh[j]) if L[j]<v<H[j]]
                if candidates:
                    _,axis,mid=max(candidates);basin_faces+=1
            b=H.copy();b[axis]=mid;a=L.copy();a[axis]=mid
            stack.extend([(a,H,depth+1),(L,b,depth+1)]);splits+=1
        live=[[*l.tolist(),*h.tolist(),n]for l,h,n in stack+stalled]
        assert sum(end.values())+len(live)==splits+1
        return dict(visits=visits,splits=splits,contractions=cuts,closed=end,material_trials=trials,basin_face_splits=basin_faces,live=live)

def main():
    if not __debug__:raise SystemExit("Assertions required")
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("checkout","cert","high","input","output","cache"):p.add_argument(key,type=Path)
    p.add_argument("--budget",type=int,default=1024);p.add_argument("--rows",help="comma-separated 0-based frontier rows")
    a=p.parse_args()
    inp=json.loads(a.input.read_bytes())
    basin=a.high/"five-body-contact-basin.json"
    assert digest(basin)=="282fafa2a4d37f72871c082b9e3e3e71c19fdb9ee4889cfa1a08fef19dcd150a"
    assert inp["target"]==prior.TARGET
    b=prior.inward_box(json.loads(basin.read_bytes()))
    source=HERE/"tomar_dynamic_affine.cpp";driver=HERE/"tomar_dynamic_five.py";mat=HERE/"tomar_five_contact_material.py"
    glob=a.cert/"SHA256SUMS";assert digest(glob)==prior.PINNED_MANIFEST
    pins={n.removeprefix("./"):v for v,n in (s.split()for s in glob.read_text().splitlines())}
    assert digest(a.cert/"part0000.json")==pins["part0000.json"]
    up=a.checkout.resolve();rev="10af98606c98ca9ffd0d9f402cae1e2efe74b792"
    for name in ("gen8.py","gen11.py"):
        assert (up/"independent"/name).read_bytes()==subprocess.check_output(["git","-C",up,"show",rev+":independent/"+name])
    sys.path.insert(0,str(up/"independent"));import gen11
    a.output.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(dir=a.output)as t:
        lib=Path(t)/"libtomar_unique.dylib"
        subprocess.run(["/usr/bin/clang++","-O3","-ftrapv","-std=c++17",f"-fmodules-cache-path={a.cache}",
            "-shared","-fPIC","-Wl,-install_name,@rpath/libtomar_unique.dylib",source,"-o",lib],check=True)
        prior.dynamic.configure(up,a.cert,lib,prior.TARGET);assert prior.dynamic.gen8.placements is gen11.placements
        solver=Step(a.high,b)
        common={"target":prior.TARGET,"input_sha256":digest(a.input),"driver_sha256":digest(Path(__file__)),
            "material_sha256":digest(mat),"localizer_sha256":digest(Path(prior.__file__)),"dynamic_sha256":digest(driver),
            "kernel_sha256":digest(lib),"basin_sha256":digest(basin),"upstream":rev,"budget":a.budget,"part":inp["part"],"pos":inp["pos"]}
        rows=inp["result"]["live"];numbers=range(len(rows)) if not a.rows else list(map(int,a.rows.split(",")))
        for i in numbers:
            output=a.output/(a.input.stem+f"-row{i}.json");assert not output.exists()
            v=rows[i];assert len(v)==23 and all(type(x)is int for x in v)
            res=solver.tree(v,a.budget);record={**common,"row":i,"root":v,"result":res}
            tmp=output.with_suffix(".partial");tmp.write_text(json.dumps(record,separators=(",",":"))+"\n");tmp.replace(output)
            print(inp["part"],inp["pos"],i,res["visits"],len(res["live"]),res["closed"],flush=True)
if __name__=="__main__":main()
