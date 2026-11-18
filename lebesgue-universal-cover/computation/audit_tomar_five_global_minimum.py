"""Verify the four complete localization layers, or recalculate them from the 501 higher-parent roots."""
import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tempfile
from zipfile import ZipFile

for var in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","VECLIB_MAXIMUM_THREADS"):
    os.environ[var]="1"
import numpy as np
from flint import arb,ctx
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE))
import locate_tomar_five_basin as initial
import locate_tomar_five_frontier as curved
import locate_tomar_five_joint_frontier as joint
from tomar_leaf_box import leaf_box

ARCHIVES={"five-basin-layers-1-3.zip":"6cc1cec098e1ca115c03b9b8c1a6619f54601c5f8573533e43ab018fa2a4094f",
          "five-basin-layer-4.zip":"550b4a800ec45f6fff8161119751ce661bba71f867971191dd57dff7f216269a"}
PARENT="d85941b19b47cda83fc86268daa50569db2d80118b7fc722d4bf90adfa457504"
PARENT_REPLAY="44390237166e1d7ee275b52a6d687d78a78f17c60236351dcbd40bb7aa866ab4"
BASIN="282fafa2a4d37f72871c082b9e3e3e71c19fdb9ee4889cfa1a08fef19dcd150a"
STATIONARY="1ad6912c1939657959fda04e6fea3312688516343ddadde9803ce70f522c377a"
WIDE="8b74545bb789b0acb8a1f17079f047942ba66fb9b6421ef56f78b587c24c6135"
UPSTREAM="10af98606c98ca9ffd0d9f402cae1e2efe74b792"
TARGET="0.83729500622102525"

def digest(file):return sha256(file.read_bytes()).hexdigest()
def patch_anchor_checks():
    # The same inner-chain area inequality uses a homotopy from the anchor.
    # Verify its endpoint inball as well as the inball at the pose under test.
    cls=curved.Material;old=cls.anchor
    def checked(self,a,b):
        data=old(self,a,b);pose=data[0]
        for j,n in enumerate((3,5,7,9)):
            assert abs(pose[3*j])+abs(pose[3*j+1])<1-1/(2*(self.st.pi/(2*n)).cos())
        return data
    cls.anchor=checked

def main():
    if not __debug__:raise SystemExit("Assertions required.")
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("checkout","cert","table","high","scratch","cache"):p.add_argument(key,type=Path)
    p.add_argument("--replay",action="store_true")
    p.add_argument("--from-layer",type=int,choices=(1,2,3,4),default=1,help="start chosen arithmetic at this layer after checking the earlier frozen coverage")
    p.add_argument("--check-moving-anchors",action="store_true",help="also verify every participating body and material primitive at all multiple contacts of each usable anchor")
    p.add_argument("--parts",help="recalculate only these root parts, while structurally auditing every layer")
    p.add_argument("--replay-parent",action="store_true")
    p.add_argument("--replay-local",action="store_true",help="recalculate stationary, transverse and both neighborhood certificates")
    a=p.parse_args();high=a.high.resolve();checkout=a.checkout.resolve()
    assert not(a.replay and a.parts)
    parent=high/"five-global-localization.json";pr=high/"five-global-localization-replay.json"
    basin=high/"five-body-contact-basin.json";stat=high/"five-body-contact-stationary.json";wide=HERE/"five-wide-contact-tube.json"
    assert digest(parent)==PARENT and digest(pr)==PARENT_REPLAY and digest(basin)==BASIN and digest(stat)==STATIONARY and digest(wide)==WIDE
    assert digest(a.table)=="fe5b816060d713642dc44e3f8981f4ae5d1b2b1c1e271b7dff75f031dd75d87c"
    par=json.loads(parent.read_bytes());rep=json.loads(pr.read_bytes());sr=json.loads(stat.read_bytes());wr=json.loads(wide.read_bytes())
    assert rep["complete"] and rep["parts"]==2048 and rep["totals"]==par["totals"]==[6710314,113881,49113935,9279,28911260,1759,32560355,501,39424759]
    assert rep["source_sha256"]==digest(high/"five_grouped_localization.cpp") and rep["driver_sha256"]==digest(high/"replay_five_global_localization.py")
    assert wr["code_sha256"]==digest(HERE/"certify_five_wide_contact_tube.py") and wr["complete"]
    assert sr["script_sha256"]==digest(high/"certify_five_body_contact_stationary.py") and arb(sr["area_ball"])<arb(TARGET)
    assert Fraction(TARGET)<Fraction(par["target"])
    m=a.cert/"SHA256SUMS";assert digest(m)==initial.PINNED_MANIFEST
    pins={n.removeprefix("./"):h for h,n in (s.split() for s in m.read_text().splitlines())}
    assert digest(a.cert/"part0000.json")==pins["part0000.json"]
    native=HERE/"tomar_dynamic_affine.cpp";dy=HERE/"tomar_dynamic_five.py";mat=HERE/"tomar_five_contact_material.py"
    joint_src=HERE/"tomar_five_joint_contact.py";st0=HERE/"locate_tomar_five_basin.py";st1=HERE/"locate_tomar_five_frontier.py";st2=HERE/"locate_tomar_five_joint_frontier.py"
    allparts=set(map(int,par["retained"]));chosen=allparts if a.replay else (set(map(int,a.parts.split(",")))if a.parts else set())
    assert chosen<=allparts
    assert subprocess.check_output(["git","-C",checkout,"rev-parse","HEAD"],text=True).strip()==UPSTREAM
    for n in ("gen8.py","gen11.py"):
        assert (checkout/"independent"/n).read_bytes()==subprocess.check_output(["git","-C",checkout,"show",UPSTREAM+":independent/"+n])
    archive={}
    for name,h in ARCHIVES.items():
        src=HERE/name;assert digest(src)==h
        with ZipFile(src)as z:
            for f in z.namelist():
                assert f not in archive;archive[f]=z.read(f)
    expected={}
    for part,ids in par["retained"].items():
        k=int(part);head=a.cert/f"part{k:04}.json";tree=head.with_suffix(".npy")
        assert digest(head)==pins[head.name]and digest(tree)==pins[tree.name]
        info=json.loads(head.read_bytes());co=np.load(tree,allow_pickle=False);pre=np.r_[1,1+np.cumsum(np.where(co<11,1,-1))]
        for pos in ids:
            c,w=leaf_box(pre,co,k,pos,info)
            root=[(Fraction(x)+s*Fraction(y))*2**44 for s in (-1,1)for x,y in zip(c,w)]
            assert all(u.denominator==1 for u in root)
            expected[f"part{k:04}-pos{pos}"]=(k,pos,list(map(int,root)),None,None)
    assert len(expected)==501
    oldbox=initial.inward_box(json.loads(basin.read_bytes()));seen=set();levels=[];count_basin=0;solver=None
    with tempfile.TemporaryDirectory(dir=a.scratch)as tmp:
        temp=Path(tmp)
        if a.replay_parent:
            subprocess.run([sys.executable,"-I","-B",high/"replay_five_global_localization.py",a.cert,temp/"parent",a.cache,"--table",a.table,"--workers","1"],check=True,stdout=subprocess.DEVNULL)
        if a.replay_local:
            for script,receipt,opts in [
                (high/"certify_five_body_contact_stationary.py",stat,[]),
                (high/"certify_five_body_contact_transverse.py",high/"five-body-contact-transverse.json",[]),
                (high/"certify_five_body_contact_basin.py",basin,[]),
                (HERE/"certify_five_wide_contact_tube.py",wide,[high,"--radius",".000012","--transverse",".00007","--cube",".00012","--cells","8","--adaptive-max","128"])]:
                # The upstream scripts import siblings from their directory.
                output=subprocess.check_output([sys.executable,"-B",script,*opts],cwd=script.parent)
                observed=json.loads(output);reference=json.loads(receipt.read_bytes())
                assert {k:v for k,v in observed.items()if k!="python"}=={k:v for k,v in reference.items()if k!="python"},script
        lib=temp/"libtomar_unique.dylib";extra=[]
        if platform.system()=="Darwin":extra.append("-Wl,-install_name,@rpath/libtomar_unique.dylib")
        compiler=shutil.which("clang++");assert compiler
        subprocess.run([compiler,"-O3","-ftrapv","-std=c++17",f"-fmodules-cache-path={a.cache}","-shared","-fPIC",*extra,native,"-o",lib],check=True)
        rebuilt_kernel=digest(lib)
        frozen_kernel="71157e1caf5583dd952931e5ae36adcbaa8403d173228e1513ecb68cb25a7cce"
        if chosen:
            sys.path.insert(0,str(checkout/"independent"));import gen11
            ctx.prec=144;initial.dynamic.configure(checkout,a.cert,lib,TARGET)
            assert initial.dynamic.gen8.placements is gen11.placements
            patch_anchor_checks()
            if a.check_moving_anchors:
                import check_five_material_anchors
                check_five_material_anchors.install()
        for layer in (1,2,3,4):
            if chosen and layer==2:solver=curved.Step(high,oldbox)
            if chosen and layer==4:solver=joint.Step(high,oldbox)
            nextset={};end=Counter();visits=splits=contracted=replayed=0;original=len(expected)
            for stem,(part,pos,roots,parent_sha,row_number)in sorted(expected.items()):
                name=f"{layer}/{stem}.json";b=archive[name];seen.add(name);r=json.loads(b)
                assert (r["part"],r["pos"])==(part,pos)and r["root"]==roots and r["target"]==TARGET
                assert r["dynamic_sha256"]==digest(dy) and r["kernel_sha256"]==frozen_kernel and r["basin_sha256"]==BASIN and r["upstream"]==UPSTREAM
                driver=(st0 if layer==1 else st1 if layer<4 else st2)
                assert r["driver_sha256"]==digest(driver)
                if layer==1:
                    assert r["source_sha256"]==digest(native) and r["tree_manifest_sha256"]==initial.PINNED_MANIFEST and r["parent_sha256"]==PARENT
                else:
                    assert r["input_sha256"]==parent_sha and r["row"]==row_number
                    assert r["localizer_sha256"]==digest(st0) and r["material_sha256"]==digest(mat)
                    if layer==4:assert r["joint_source_sha256"]==digest(joint_src) and r["wide_tube_sha256"]==WIDE
                budget=[4096,768,1536,1536][layer-1];assert r["budget"]==budget
                result=r["result"];live=result["live"];ends=result["closed"]
                allowed={"radial","basin","area1","area2","area3"}|({"material"}if layer>1 else set())
                assert set(ends)<=allowed and all(type(v)is int and v>0 for v in ends.values())
                s=result["splits"];c=result["contractions"];v=result["visits"]
                assert 1<=v<=budget and s>=0 and c>=0 and sum(ends.values())+len(live)==s+1
                stall=v-s-c-sum(ends.values());assert 0<=stall<=len(live) and (v==budget or stall==len(live))
                for j,cell in enumerate(live):
                    assert len(cell)==23 and all(type(u)is int for u in cell) and cell[-1]>=0
                    assert all(roots[i]<=cell[i]<=cell[11+i]<=roots[11+i]for i in range(11))
                    nextset[f"{stem}-row{j}"]=(part,pos,cell,sha256(b).hexdigest(),j)
                if part in chosen and layer>=a.from_layer:
                    actual=initial.cover(roots,oldbox,budget)if layer==1 else solver.tree(roots,budget)
                    assert actual==result,(layer,stem)
                    replayed+=1
                end.update(ends);visits+=v;splits+=s;contracted+=c
            levels.append(dict(layer=layer,input_boxes=original,visits=visits,splits=splits,contractions=contracted,
                terminals=dict(sorted(end.items())),remaining=len(nextset),original_roots_remaining=sorted(set((v[0],v[1])for v in nextset.values())),replayed=replayed))
            count_basin+=end["basin"];expected=nextset
    assert not expected and set(archive)==seen and count_basin==10
    ctx.prec=128
    base=json.loads((a.cert/"part0000.json").read_bytes());outside=[]
    for i,n in enumerate(base["ns"]):
        k=3*i;inner=1-1/(2*(arb.pi()/(2*n)).cos())
        d=min(-base["lo"][k],base["hi"][k],-base["lo"][k+1],base["hi"][k+1])
        v=arb.pi()*(arb(1)/4+inner*inner)/2+arb(d)*(arb(1)/2+inner)
        assert v>arb(par["target"]);outside.append(str(v))
    print(json.dumps(dict(complete=True,finite_family="disk and regular Reuleaux 3,5,7,9-gons",unique_modulo_rigid_symmetries=True,
        target=TARGET,minimum_interval=sr["area_ball"],higher_parent=PARENT,higher_parent_complete_replay=PARENT_REPLAY,
        original_roots=2048,canonical_leaves=6710314,higher_target_exclusions=6709813,exceptions=501,basin_terminals=count_basin,
        levels=levels,archives=ARCHIVES,replayed_parts=sorted(chosen),parent_replayed=a.replay_parent,local_certificates_replayed=a.replay_local,
        full_moving_anchor_classification=a.check_moving_anchors,replayed_start_layer=a.from_layer,rebuilt_kernel_sha256=rebuilt_kernel,kernel_source_sha256=digest(native),root_exterior=outside,historical_area_inputs=[]),indent=2))

if __name__=="__main__":main()
