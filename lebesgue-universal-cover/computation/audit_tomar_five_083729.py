"""Verify or replay the strict five-body 0.83729 bound on every exception of the high-target parent."""
import argparse
from fractions import Fraction as F
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from zipfile import ZipFile
for k in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","VECLIB_MAXIMUM_THREADS"):os.environ[k]="1"
import numpy as np
from flint import arb,ctx
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from tomar_d3_mask import leaf_mask
from tomar_leaf_box import leaf_box
from verify_tomar_fan_table import PINNED_MANIFEST
import tomar_dynamic_five as dynamic

def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def digest(b):return hashlib.sha256(b).hexdigest()

def main():
    if not __debug__:raise SystemExit("Verification requires enabled assertions.")
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("checkout","cert","table","high_dir","archive","scratch","cache"):p.add_argument(k,type=Path)
    p.add_argument("--replay-parent",action="store_true")
    p.add_argument("--replay",action="store_true")
    p.add_argument("--replay-parts",help="recompute just these high-exception roots in addition to the complete structural audit")
    a=p.parse_args();high=a.high_dir.resolve();target="0.83729"
    assert sha(a.archive)=="d040caaf37c8d6dad2bc2273f5cc6364a6862b9f5696b6de666f8276a65ab6b7"
    hfile=high/"five-global-localization.json";proof=high/"five-global-localization-replay.json";ixfile=high/"five-orbit-slice.json"
    assert sha(hfile)=="d85941b19b47cda83fc86268daa50569db2d80118b7fc722d4bf90adfa457504" and sha(proof)=="44390237166e1d7ee275b52a6d687d78a78f17c60236351dcbd40bb7aa866ab4"
    assert sha(ixfile)=="ec483a5d64c86fc1d305c63081e08d66fd128bb2033622764ee26b2d45377e02" and sha(a.table)=="fe5b816060d713642dc44e3f8981f4ae5d1b2b1c1e271b7dff75f031dd75d87c"
    parent=json.loads(hfile.read_bytes());hp=json.loads(proof.read_bytes());ix=json.loads(ixfile.read_bytes())
    assert hp["complete"] and hp["parts"]==2048 and hp["totals"]==parent["totals"]==[6710314,113881,49113935,9279,28911260,1759,32560355,501,39424759]
    assert hp["source_sha256"]==sha(high/"five_grouped_localization.cpp") and hp["driver_sha256"]==sha(high/"replay_five_global_localization.py")
    full={(int(k),i)for k,v in parent["retained"].items()for i in v};assert len(full)==501 and F(target)<F(parent["target"])
    all_exception_roots={k for k,i in full}
    assert not (a.replay and a.replay_parts)
    replay_roots=all_exception_roots if a.replay else set(map(int,a.replay_parts.split(","))) if a.replay_parts else set()
    assert replay_roots<=all_exception_roots
    raw_manifest=a.cert/"SHA256SUMS";assert sha(raw_manifest)==PINNED_MANIFEST
    pinned={n.removeprefix("./"):h for h,n in (s.split()for s in raw_manifest.read_text().splitlines())}
    rev="10af98606c98ca9ffd0d9f402cae1e2efe74b792";checkout=a.checkout.resolve()
    assert subprocess.check_output(["git","-C",checkout,"rev-parse","HEAD"],text=True).strip()==rev
    for name in ("gen8.py","gen11.py"):
        x=checkout/"independent"/name
        assert digest(subprocess.check_output(["git","-C",checkout,"show",rev+":independent/"+name]))==sha(x)
    sys.path.insert(0,str(checkout/"independent"))
    import gen8
    default_placements=gen8.placements
    import gen11
    five_placements=gen11.placements
    assert len(default_placements([0]*11))==3 and len(five_placements([0]*11))==4
    source=HERE/"tomar_dynamic_affine.cpp"
    with tempfile.TemporaryDirectory(dir=a.scratch)as tmp,ZipFile(a.archive)as z:
        d=Path(tmp)
        if a.replay_parent:
            subprocess.run([sys.executable,"-I","-B",high/"replay_five_global_localization.py",a.cert,d/"high",a.cache,
                "--table",a.table,"--workers","1"],stdout=subprocess.DEVNULL,check=True)
        lib=d/"libtomar_five_083729.dylib"
        subprocess.run(["/usr/bin/clang++","-O3","-ftrapv","-std=c++17",f"-fmodules-cache-path={a.cache}",
            "-shared","-fPIC","-Wl,-install_name,@rpath/libtomar_five_083729.dylib",source,"-o",lib],check=True)
        dynamic.configure(checkout,a.cert,lib,target)
        active_raw=z.read("active.json");act=json.loads(active_raw);assert act["high_selection_sha256"]==sha(hfile)
        assert act["target"]==target and act["parent_target"]==parent["target"]
        first={(part,pos):h for part,pos,h in act["first_receipts"]}
        assert len(first)==len(act["first_receipts"]) and set(first)<=full
        nextset={(int(k),pos)for k,ids in act["retained"].items()for pos in ids}
        assert set().union(*(set(group)for group in act["assignment"]))==set(map(int,act["retained"]))
        assert sum(map(len,act["assignment"]))==len(act["retained"])
        seen={"active.json"};accepted=set();counts={};cache={};nreplayed=0;initial_failed=[]
        for stage,ids in (("first",set(first)),("five",nextset)):
            dynamic.gen8.placements=default_placements if stage=="first" else five_placements
            if stage=="five":assert nextset==full-accepted and len(accepted)==act["first_closed"]
            visits=accepted_here=total=0
            for part,pos in sorted(ids):
                name=f"{stage}/part{part:04}-pos{pos}.json";seen.add(name);data=z.read(name);record=json.loads(data)
                if stage=="first":assert digest(data)==first[(part,pos)]
                if part not in cache:
                    stem=f"part{part:04}";hf=a.cert/(stem+".json");tr=a.cert/(stem+".npy")
                    assert sha(hf)==pinned[hf.name] and sha(tr)==pinned[tr.name]
                    header=json.loads(hf.read_bytes());codes=np.load(tr,allow_pickle=False);cls=leaf_mask(codes,part,F(header["hi"][2]))
                    assert hashlib.sha256(cls.tobytes()).hexdigest()==ix["rows"][part][-1] and header["part"]==part
                    pref=np.r_[1,1+np.cumsum(np.where(codes<11,1,-1))]
                    cache[part]=header,codes,cls,pref,pinned[hf.name],pinned[tr.name]
                header,codes,cls,pref,hhead,htree=cache[part]
                assert cls[pos] in (1,2)
                drv=HERE/("run_tomar_dynamic_083729.py"if stage=="first"else"run_tomar_five_proposals_083729.py")
                expected={"target":target,"selection_sha256":sha(hfile),"slice_index_sha256":sha(ixfile),"input_manifest_sha256":PINNED_MANIFEST,
                    "source_sha256":sha(source),"library_sha256":sha(lib),"affine_driver_sha256":sha(HERE/"tomar_dynamic_five.py"),
                    "driver_sha256":sha(drv),"decoder_sha256":sha(HERE/"tomar_leaf_box.py"),"proposal_sha256":sha(checkout/"independent/gen8.py"),
                    "upstream":rev,"moves_seed":20261118,"max_calls":8192,"max_depth":96,"moves":65,
                    "part":part,"header_sha256":hhead,"codes_sha256":htree}
                if stage=="five":
                    expected.update(active_sha256=digest(active_raw),five_body_proposal_sha256=sha(checkout/"independent/gen11.py"))
                assert all(record[k]==v for k,v in expected.items())
                q=record["row"];assert q["id"]==part and q["label"]==pos and 1<=q["vis"]<=8192 and 0<=q["md"]<=96
                total+=1
                assert set(q["closed"])<={"1","2","3"} and all(type(v)is int and v>0 for v in q["closed"].values())
                assert q["splits"]>=0 and q["contractions"]>=0
                if q["remaining"]:
                    assert stage=="first" and (part,pos)in nextset
                    initial_failed.append([part,pos,len(q["remaining"])])
                    continue
                term=sum(q["closed"].values());assert term==q["splits"]+1 and q["vis"]==q["splits"]+q["contractions"]+term
                assert (part,pos)not in accepted;accepted.add((part,pos));accepted_here+=1;visits+=q["vis"]
                if part in replay_roots:
                    c,w=leaf_box(pref,codes,part,pos,header)
                    v=[(F(x)-F(y))*2**44 for x,y in zip(c,w)]+[(F(x)+F(y))*2**44 for x,y in zip(c,w)]
                    assert all(x.denominator==1 for x in v)
                    actual=dynamic.work((part,pos,list(map(int,v)),8192,96))
                    assert json.loads(json.dumps(actual))==q
                    nreplayed+=1
            counts[stage]={"tested":total,"accepted":accepted_here,"successful_box_visits":visits}
        assert sorted(initial_failed)==sorted(act["first_incomplete"])
        assert accepted==full and len(z.namelist())==len(seen) and set(z.namelist())==seen
        assert nreplayed==sum(k in replay_roots for k,i in full)
    D=2**57;v=F(target)*D;N=-(-v.numerator//v.denominator)
    ctx.prec=128;root_file=a.cert/"part0000.json"
    assert sha(root_file)==pinned[root_file.name]
    base=json.loads(root_file.read_bytes());outside=[]
    for i,n in enumerate(base["ns"]):
        k=3*i;r=1-1/(2*(arb.pi()/(2*n)).cos());dist=min(-base["lo"][k],base["hi"][k],-base["lo"][k+1],base["hi"][k+1])
        val=arb.pi()*(arb(1)/4+r*r)/2+arb(dist)*(arb(1)/2+r);assert val>arb(target);outside.append(str(val))
    print(json.dumps({"target":target,"complete":True,"canonical_leaves":6710314,"higher_threshold_closed":6709813,
        "original_exceptions":501,"remaining":0,"area_numerator":N,"area_denominator":D,"levels":counts,
        "successful_original_trees_replayed":nreplayed,"replayed_exception_root_parts":sorted(replay_roots),"higher_parent_replayed_here":a.replay_parent,"root_exterior":outside,
        "higher_parent_sha256":sha(hfile),"higher_parent_complete_replay_sha256":sha(proof),"historical_area_inputs":[],
        "archive_sha256":sha(a.archive),"archive_members":len(seen),"archive_bytes":a.archive.stat().st_size},indent=2))

if __name__=="__main__":main()
