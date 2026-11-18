"""Replay the .83744 weighted point trees and their complementary five-body certificates."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import ctypes
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import platform
import shutil
import subprocess
import sys
from zipfile import ZipFile

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import flint
from flint import arb, ctx
import numpy as np
from certify_nine import ar
from certify_tomar_six_high import FINE, Q, SourceTree, build_kernel, ceil, load_table
from certify_finite_point_leaf import finite_clouds
from replay_finite_point_083732 import PointProof

TARGET = F(83744, 100000)


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def read(path, expected):
    assert digest(path)==expected
    return json.loads(path.read_bytes())


def scaled(values):
    z = [F(v)*FINE for v in values]
    assert all(q.denominator==1 for q in z)
    return np.array(list(map(int,z)),dtype=np.int64)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("tree","table32","table64","work","cache","retained","outer_audit"):
        p.add_argument(name,type=Path)
    p.add_argument("--workers",type=int,default=3)
    p.add_argument("--five-kit",type=Path)
    p.add_argument("--upstream",type=Path)
    p.add_argument("--dense-source",type=Path,default=HERE)
    p.add_argument("--parts",nargs="+",type=int)
    args=p.parse_args()
    assert __debug__ and 1<=args.workers<=32
    assert (args.five_kit is None)==(args.upstream is None)
    work=args.work.resolve()
    assert not work.exists()
    rec_path=HERE/"finite-point-lift-083744.json"
    rec=json.loads(rec_path.read_bytes())
    roots=read(args.retained, rec["retained_sha256"])
    audit=read(args.outer_audit, rec["outer_audit_sha256"])
    assert F(rec["target"])==TARGET==F(audit["target"]) and audit["complete"]
    goal=ceil(2*Q*Q*TARGET)
    assert audit["integer_goal"]==goal==120687823102244742
    assert audit["parts"]==1024 and len(audit["audits"])==4
    all_keys={(int(part),pos) for part,ids in roots.items() for pos in ids}
    old=read(HERE/"tomar-083743-retained.json",rec["previous_sha256"])
    prior={(int(p),i) for p,ids in old.items() for i in ids}
    assert len(prior)==2659 and len(prior-all_keys)==rec["dropped"] and len(all_keys-prior)==rec["added"]
    assert len(all_keys)==rec["retained_count"]==audit["audits"][-1]["remaining"]
    skipped={tuple(r) for r in rec["five_only"]}
    assert skipped<=all_keys and len(skipped)==len(rec["five_only"])
    dense=rec["dense"]
    dense_keys=set(map(tuple,dense["original_leaves"]))
    assert len(dense_keys)==len(dense["original_leaves"]) and dense_keys<=all_keys-skipped and (180,24225) in dense_keys
    ordinary_keys=all_keys-skipped-dense_keys
    rows=rec["rows"]
    assert [tuple(r[:2]) for r in rows]==sorted(ordinary_keys)

    cloud=read(HERE/"finite-point-leaf-input.json",rec["point_sha256"])
    bank=read(HERE/"finite-point-global-input.json",rec["bank_sha256"])
    base=read(HERE/"tomar-six-high-input.json",rec["base_sha256"])
    ctx.prec=160
    r=ar(F(707,1000))
    assert r*r<ar(F(1,2)) and (arb.pi()-(1/(2*r)).acos())/4+(r*r-ar(F(1,4))).sqrt()/2>ar(F(".838"))
    tree=SourceTree(args.tree)
    tables={n:load_table(getattr(args,f"table{n}"),tree,base["tables"][str(n)]) for n in (32,64)}
    point=PointProof(cloud,bank["cells"],100000000008865717121)
    cache=args.cache.resolve()
    _, compiler=build_kernel(work,base,cache)
    src=HERE/"tomar_six_083730_guard.cpp"
    assert digest(src)==rec["weighted_sha256"]=="42a4ee583d22940b3cd0b05f04285ba7e3d351745de707476211dce4efba5f81"
    library=work/"guard.so"
    clang=shutil.which("clang++")
    assert clang
    subprocess.run([clang,"-O3","-std=c++17","-ftrapv","-shared","-fPIC",
                    f"-fmodules-cache-path={cache}","-I",str(work),str(src),"-o",str(library)],check=True)
    guard=ctypes.CDLL(str(library)).lift_guard_run
    v=lambda d:np.ctypeslib.ndpointer(dtype=d,ndim=1,flags="C_CONTIGUOUS")
    i64,i32=v(np.int64),v(np.int32)
    guard.argtypes=[i64,i64,i64,i32,i32,ctypes.c_int64,ctypes.c_int,ctypes.c_int,
                    i64,i64,i64,i64,i64,ctypes.c_int,ctypes.c_int,ctypes.c_int,i64]
    guard.restype=ctypes.c_int
    root=scaled(tree.root["hi"])

    def verify(row, box, support):
        part,pos,n,calls,old_count,new_count,depth,tried,joint=row
        lo,hi=(scaled(z) for z in zip(*box,strict=True))
        st=np.zeros(30,dtype=np.int64)
        assert guard(lo,hi,root,*tables[n],goal,depth+1,calls,*support,
                     len(support[1]),len(support[3])-1,2,st)==0
        assert st[:7].tolist()==[1,calls,old_count,new_count,depth,0,tried] and st[29]==joint,(part,pos,st[:7])
        assert calls==2*(old_count+new_count)-1
        return row

    selected=[r for r in rows if args.parts is None or r[0] in args.parts]
    # Table size scales hull work; wide failed boxes also evaluate the cap inequalities.
    # The score orders tasks only. Every node decision is checked against its full record.
    def weight(row):
        return row[3]*sum(int(n) for n in tables[row[2]][1]) + row[7]*105
    ordered=sorted(selected,key=weight,reverse=True)
    with ThreadPoolExecutor(args.workers) as ex:
        pending={}
        for r in ordered:
            box=tree.leaf(*r[:2])
            pending[ex.submit(verify,r,box,point.presolve(box))]=tuple(r[:2])
        results=sorted((f.result() for f in as_completed(pending)),key=lambda r:tuple(r[:2]))
    assert results==selected
    out=dict(target=str(TARGET),ordinary_complete=len(results)==len(ordinary_keys),terminal_trees=len(results),ranked_work_units=sum(weight(r) for r in results),
             nodes=sum(r[3] for r in results),five_terminals=sum(r[4] for r in results),
             point_terminals=sum(r[5] for r in results),joint_terminals=sum(r[8] for r in results),
             point_used=sum(r[5]>0 for r in results),pure_five_replayed=False,complete=False,
             ordinary_input_sha256=digest(rec_path),retained_sha256=digest(args.retained),
             outer_audit_sha256=digest(args.outer_audit),compiler=compiler,python=platform.python_version(),
             flint=flint.__version__,driver_sha256=digest(Path(__file__)))
    if args.five_kit:
        assert args.parts is None and digest(args.five_kit)==rec["five_kit_sha256"]
        kit=work/"five-kit"
        with ZipFile(args.five_kit) as z:
            assert all(not Path(n).is_absolute() and ".." not in Path(n).parts for n in z.namelist())
            z.extractall(kit)
        sel=json.loads((kit/"expected/selection.json").read_bytes())
        assert {(r[0],r[1]) for r in sel["leaves"]}==skipped and len(sel["leaves"])==len(skipped)
        run=work/"five-replay"
        subprocess.run([sys.executable,"-I","-B",kit/"replay_083744.py",args.upstream.resolve(),
                        args.tree.resolve(),run,cache,"--workers",str(min(args.workers,12))],check=True)
        pure=[]
        for part,pos in sorted(skipped):
            stem=f"part{part}-pos{pos}"
            output=run/stem/(stem+".json")
            d=json.loads(output.read_bytes())
            assert (d["part"],d["pos"])==(part,pos) and F(d["threshold"])>=TARGET and not d["result"]["remaining"]
            pure.append([part,pos,d["result"]["vis"],d["trace_sha256"]])
        dense_comp=args.dense_source.resolve()
        shared_comp=Path(sys.modules["certify_tomar_six_high"].__file__).resolve().parent
        dense_file=dense_comp/dense["filename"]
        assert digest(dense_file)==dense["sha256"]
        dense_data=json.loads(dense_file.read_bytes())
        assert dense_data["original_leaves"]==dense["original_leaves"]
        assert digest(dense_comp/"check_scaled_p30_leaves_083744.py")==dense["driver_sha256"]
        dense_out=json.loads(subprocess.check_output([
            sys.executable,"-I","-B",dense_comp/"check_scaled_p30_leaves_083744.py",shared_comp,args.tree,
            dense_file,"--replay",args.table64,work/"dense-replay",cache,"--workers",str(args.workers)]))
        assert dense_out["complete"] and dense_out["geometric_replayed"]
        assert dense_out["leaves"]==len(dense_keys) and dense_out["visits"]==dense["visits"]
        out["dense_replayed"]=dense_out
        out.update(pure_five_replayed=True,pure_leaves=pure,kit_sha256=rec["five_kit_sha256"],complete=out["ordinary_complete"])
    if out["complete"]:
        groups,e,_=finite_clouds()
        groups.append(dict(body="P",denominator=cloud["denominator"],points=cloud["points"],
                           diameter_squared_numerator=100000000008865717121))
        path=work/"six-rational-groups.json"
        path.write_text(json.dumps(groups,separators=(",",":"))+"\n")
        assert digest(path)=="6a264f46e692045b6fd7b5c74f33f86b71601e972f9236479962a5d26742be4b"
        floor=(TARGET-F(22,7)*e**2)/(1+2*e/(F(1,2)-e))
        assert floor>F(".83743015029")
        out.update(all_finite_area_lower=str(floor),all_finite_decimal_floor=".83743015029",cloud_sha256=digest(path))
    (work/"replay.json").write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps(out,indent=2),flush=True)


if __name__=="__main__":
    main()
