"""Verify the high-threshold parent and a complete cover of its 501 exceptional original leaves."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from zipfile import ZipFile
import numpy as np
from flint import arb,ctx

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from tomar_d3_mask import leaf_mask
from verify_tomar_fan_table import PINNED_MANIFEST,verified_table

def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def run(s,*argv):return subprocess.check_output([sys.executable,"-I","-B",s,*map(str,argv)])

p=argparse.ArgumentParser(description=__doc__)
for k in ("checkout","cert","table","high_dir","archive","scratch","cache"):p.add_argument(k,type=Path)
p.add_argument("--replay-parent",action="store_true")
p.add_argument("--replay-stages",action="store_true")
p.add_argument("--replay-dynamic",action="store_true")
a=p.parse_args();high=a.high_dir.resolve();target="0.83725"
highfile=high/"five-global-localization.json";proof=high/"five-global-localization-replay.json";index=high/"five-orbit-slice.json"
assert sha(highfile)=="d85941b19b47cda83fc86268daa50569db2d80118b7fc722d4bf90adfa457504"
assert sha(proof)=="44390237166e1d7ee275b52a6d687d78a78f17c60236351dcbd40bb7aa866ab4"
assert sha(index)=="ec483a5d64c86fc1d305c63081e08d66fd128bb2033622764ee26b2d45377e02"
assert sha(a.table)=="fe5b816060d713642dc44e3f8981f4ae5d1b2b1c1e271b7dff75f031dd75d87c"
sp=high/"replay_five_global_localization.py";first_kernel=high/"five_grouped_localization.cpp"
hd=json.loads(highfile.read_bytes());hp=json.loads(proof.read_bytes());ix=json.loads(index.read_bytes())
assert hp["complete"] and hp["parts"]==2048 and hp["totals"]==hd["totals"]==[6710314,113881,49113935,9279,28911260,1759,32560355,501,39424759]
assert hp["source_sha256"]==sha(first_kernel) and hp["driver_sha256"]==sha(sp)
assert F(hd["target"])==F("0.837295006222")>F(target)
parent_left={(int(k),i)for k,pos in hd["retained"].items()for i in pos}
assert len(parent_left)==501 and len(hd["retained"])==19
manifest=a.cert/"SHA256SUMS";assert sha(manifest)==PINNED_MANIFEST
pinned={n.removeprefix("./"):h for h,n in (l.split()for l in manifest.read_text().splitlines())}
base=json.loads((a.cert/"part0000.json").read_bytes());verified_table(a.table,a.cert,base)
ctx.prec=128;exterior=[]
for i,n in enumerate(base["ns"]):
    k=3*i;r=1-1/(2*(arb.pi()/(2*n)).cos());dist=min(-base["lo"][k],base["hi"][k],-base["lo"][k+1],base["hi"][k+1])
    ar=arb.pi()*(arb(1)/4+r*r)/2+arb(dist)*(arb(1)/2+r)
    assert ar>arb(target);exterior.append(str(ar))
fine=HERE/"tomar_fine_angle_fan.cpp";dyn=HERE/"tomar_dynamic_affine.cpp"
limits=[(65535,55),(1048575,60),(1048575,100),(8388607,115)]
known_selections=[1759,221,129,22,16];levels=[];before=None;prior_dir=None
with tempfile.TemporaryDirectory(dir=a.scratch)as folder,ZipFile(a.archive)as archive:
    temp=Path(folder);names=set()
    if a.replay_parent:
        run(sp,a.cert,temp/"high-replay",a.cache,"--table",a.table,"--workers",1)
    binaries={}
    for source,name in [(first_kernel,"libtomar_localized_target.dylib"),(fine,"libtomar_localized_target.dylib"),(dyn,"libtomar_dynamic_five.dylib")]:
        out=temp/name
        subprocess.run(["/usr/bin/clang++","-O3","-ftrapv","-std=c++17",f"-fmodules-cache-path={a.cache}","-shared","-fPIC",
            "-Wl,-install_name,@rpath/"+name,source,"-o",out],check=True);binaries[sha(source)]=sha(out)
    for j in range(1,6):
        prefix=f"level{j}/";sd=temp/f"level{j}";sd.mkdir();rf=sd/"selection.json"
        names.add(prefix+rf.name);rf.write_bytes(archive.read(prefix+rf.name));s=json.loads(rf.read_bytes())
        assert sum(map(len,s["retained"].values()))==known_selections[j-1]
        if j==1:
            assert sha(rf)=="023efd50c4aa3ac465339ecae7998c4e79a97a2884c7160e4d4c9ac993fe799f"
            assert parent_left<={(int(k),i)for k,pos in s["retained"].items()for i in pos}
        else:
            regenerated=temp/f"selection-{j}.json"
            run(HERE/"select_tomar_localized_residues.py",before,prior_dir,regenerated)
            assert rf.read_bytes()==regenerated.read_bytes()
        rc=temp/f"repeat-{j}"
        kernel=first_kernel if j<3 else fine
        driver=HERE/("check_tomar_localized_target.py" if j<5 else "replay_tomar_dynamic_exceptions.py")
        if j<5 and a.replay_stages:
            n,h=limits[j-1]
            run(driver,rf,a.cert,a.table,index,kernel,rc,a.cache,"--target",target,"--calls",n,"--depth",h)
        if j==5 and a.replay_dynamic:
            run(driver,a.checkout,a.cert,index,rf,rc,a.cache)
        accounted=set();checked=calls=0
        for k,pos in sorted(s["retained"].items(),key=lambda kv:int(kv[0])):
            part=int(k);filename=f"part{part:04}.json";f=sd/filename;n=prefix+filename
            names.add(n);f.write_bytes(archive.read(n));r=json.loads(f.read_bytes())
            hf=a.cert/filename;cf=hf.with_suffix(".npy")
            assert sha(hf)==pinned[hf.name] and sha(cf)==pinned[cf.name]
            h=json.loads(hf.read_bytes());codes=np.load(cf,allow_pickle=False)
            classes=leaf_mask(codes,part,F(h["hi"][2]));assert hashlib.sha256(classes.tobytes()).hexdigest()==ix["rows"][part][-1]
            assert pos==sorted(set(pos)) and np.isin(classes[pos],(1,2)).all()
            for key,expected in {"part":part,"target":target,"selection_sha256":sha(rf),"slice_index_sha256":sha(index),
                    "header_sha256":sha(hf),"codes_sha256":sha(cf),"driver_sha256":sha(driver),
                    "source_sha256":sha(kernel if j<5 else dyn),"library_sha256":binaries[sha(kernel if j<5 else dyn)]}.items():
                assert r[key]==expected,(j,part,key)
            if j<5:
                n,h=limits[j-1]
                assert r["depth"]==h and r["calls"]==n and r["table_sha256"]==sha(a.table)
                rest=r["remaining"];assert rest==sorted(set(rest)) and set(rest)<=set(pos)
                assert r["tested"]==len(pos) and r["closed"]==len(pos)-len(rest) and r["maximum_depth"]<=h
                accounted.update((part,i)for i in set(pos)-set(rest));calls+=r["checked_boxes"]
                if a.replay_stages:assert (rc/filename).read_bytes()==f.read_bytes()
            else:
                assert r["affine_driver_sha256"]==sha(HERE/"tomar_dynamic_five.py") and r["upstream"]=="10af98606c98ca9ffd0d9f402cae1e2efe74b792"
                assert r["proposal_sha256"]==sha(a.checkout/"independent/gen8.py") and r["decoder_sha256"]==sha(HERE/"tomar_leaf_box.py")
                assert r["max_calls"]==8192 and r["max_depth"]==90 and r["moves"]==65 and r["moves_seed"]==20261118
                assert [z["label"]for z in r["rows"]]==pos
                for z in r["rows"]:
                    assert z["id"]==part and not z["remaining"] and z["md"]<=90 and z["vis"]<=8192
                    terminal=sum(z["closed"].values());assert terminal==z["splits"]+1
                    assert z["vis"]==terminal+z["splits"]+z["contractions"]
                    calls+=z["vis"];accounted.add((part,z["label"]))
                if a.replay_dynamic:assert (rc/filename).read_bytes()==f.read_bytes()
            checked+=len(pos)
        prev=len(parent_left);parent_left-=accounted
        levels.append({"level":j,"original_positions_tested":checked,"boxes":calls,"high_exceptions_closed":prev-len(parent_left),"high_remaining":len(parent_left)})
        before,prior_dir=rf,sd
    assert len(archive.namelist())==len(names) and set(archive.namelist())==names
assert not parent_left and [g["high_exceptions_closed"]for g in levels]==[280,92,107,6,16]
D=2**57;threshold=F(target)*D;N=-(-threshold.numerator//threshold.denominator)
print(json.dumps({"target":target,"complete":True,"canonical_leaves":6710314,"full_original_leaves":16193402,
    "prior_high_closed":6710314-501,"prior_high_target":hd["target"],"remaining":0,"levels":levels,
    "root_exterior":exterior,"integer_numerator":N,"integer_denominator":D,"historical_area_inputs":[],
    "high_input_sha256":sha(highfile),"high_complete_replay_sha256":sha(proof),"high_replayed_here":a.replay_parent,
    "fixed_stages_replayed_here":a.replay_stages,"dynamic_replayed_here":a.replay_dynamic,
    "archive_sha256":sha(a.archive),"archive_bytes":a.archive.stat().st_size,"archive_members":len(names)},indent=2))
