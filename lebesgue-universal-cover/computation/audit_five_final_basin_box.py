"""Replace the only fourth-layer basin terminal by a one-node strict material-area exclusion."""
import argparse
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tempfile
from zipfile import ZipFile

for k in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","VECLIB_MAXIMUM_THREADS"):os.environ[k]="1"
import numpy as np
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE))
from audit_tomar_five_global_minimum import ARCHIVES,BASIN,UPSTREAM,digest,patch_anchor_checks
import locate_tomar_five_joint_frontier as joint

def main():
    if not __debug__:raise SystemExit("Assertions required.")
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("checkout","cert","high","scratch","cache"):p.add_argument(key,type=Path)
    a=p.parse_args();up=a.checkout.resolve();q=a.high.resolve()
    name="five-basin-layer-4.zip";archive=HERE/name
    assert digest(archive)==ARCHIVES[name]
    sole=[]
    with ZipFile(archive)as z:
        for n in z.namelist():
            x=json.loads(z.read(n));count=x["result"]["closed"].get("basin",0)
            if count:sole.append((n,x,count))
    assert len(sole)==1
    entry,source,count=sole[0]
    assert count==1 and entry=="4/part0180-pos23955-row7-row6-row4.json"
    assert digest(q/"five-body-contact-basin.json")==BASIN
    m=a.cert/"SHA256SUMS";assert digest(m)==joint.prior.PINNED_MANIFEST
    pins={n.removeprefix("./"):v for v,n in (r.split()for r in m.read_text().splitlines())}
    assert digest(a.cert/"part0000.json")==pins["part0000.json"]
    for n in ("gen8.py","gen11.py"):
        assert (up/"independent"/n).read_bytes()==subprocess.check_output(["git","-C",up,"show",UPSTREAM+":independent/"+n])
    sys.path.insert(0,str(up/"independent"));import gen11
    with tempfile.TemporaryDirectory(dir=a.scratch)as t:
        lib=Path(t)/"libtomar_unique.dylib";flags=["-shared","-fPIC"]
        if platform.system()=="Darwin":flags+=["-Wl,-install_name,@rpath/libtomar_unique.dylib"]
        compiler=shutil.which("clang++");assert compiler
        subprocess.run([compiler,"-O3","-ftrapv","-std=c++17",f"-fmodules-cache-path={a.cache}",*flags,HERE/"tomar_dynamic_affine.cpp","-o",lib],check=True)
        joint.prior.dynamic.configure(up,a.cert,lib,joint.prior.TARGET);assert joint.prior.dynamic.gen8.placements is gen11.placements
        patch_anchor_checks()
        solver=joint.Step(q,joint.prior.inward_box(json.loads((q/"five-body-contact-basin.json").read_bytes())))
        terminal=[];old_decision=solver.inside
        def trace(l,h):
            accept=old_decision(l,h)
            if accept:terminal.append([*l.tolist(),*h.tolist(),0])
            return accept
        solver.inside=trace
        assert solver.tree(source["root"],source["budget"])==source["result"] and len(terminal)==1
        # The exact-area check must succeed with every possible basin test disabled.
        solver.inside=lambda l,h:False
        check=solver.tree(terminal[0],1)
        assert check["visits"]==1 and check["closed"]=={"material":1} and not check["live"]
    print(json.dumps(dict(stage4_archive_sha256=ARCHIVES[name],source=entry,stage4_basin_terminals=1,
        terminal=terminal[0][:22],extra_result=check,additional_basin_assumptions=[]),indent=2))
if __name__=="__main__":main()
