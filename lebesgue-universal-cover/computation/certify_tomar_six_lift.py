"""Replay the final 501+40 original leaves for the strict .83730 six-body bound."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import ctypes
from fractions import Fraction as F
import json
import platform
from pathlib import Path
import shutil
import subprocess

from flint import arb,ctx
import flint
import numpy as np
from certify_tomar_six_high import (
    HERE, FINE, Q, SourceTree, DirectionProof, build_kernel, digest, ceil, load_table
)
from certify_nine import ar


def selected(rows):
    return {(int(p),i) for p,ids in rows.items() for i in ids}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("tree","localization","outer40","table32","table64","work"):
        p.add_argument(k,type=Path)
    p.add_argument("--clang-cache",type=Path)
    p.add_argument("--workers",type=int,default=3)
    p.add_argument("--compare-caps",action="store_true")
    p.add_argument("--parts",nargs="+",type=int)
    a=p.parse_args()
    older=(HERE/"tomar-six-high-input.json").read_bytes()
    source=(HERE/"tomar-six-083730-input.json").read_bytes()
    spec=json.loads(source)
    assert digest(older)==spec["previous_proof_sha256"]
    previous=json.loads(older)
    loc=json.loads(a.localization.read_bytes())
    out=json.loads(a.outer40.read_bytes())
    assert F(out["target"])==F(spec["target"])==F(".83730")
    assert F(loc["target"])==F(previous["target"]) and loc["denominator"]==FINE
    assert digest(json.dumps(loc["retained"],sort_keys=True,separators=(",",":")).encode())==previous["localization_retained_sha256"]
    assert out["retained"]==spec["outside40"]
    K,B=selected(loc["retained"]),selected(out["retained"])
    assert len(K)==501 and len(B)==40 and not K&B
    assert [(p,i) for p,i,*_ in spec["cases"]]==sorted(K|B)
    raw=(HERE/"large-support-harmonic-input.json").read_bytes()
    assert digest(raw)==previous["design_sha256"]
    ctx.prec=160
    r=ar(F(707,1000)); pi=arb.pi()
    assert r*r<ar(F(1,2))
    assert (pi-(1/(2*r)).acos())/4+(r*r-ar(F(1,4))).sqrt()/2>ar(F(".838"))
    assert F(spec["target"])<F(".838")
    scratch=a.work.resolve()
    cache=a.clang_cache or scratch/"clang-cache"
    _,clang=build_kernel(scratch,previous,cache)
    compiler=shutil.which("clang++")
    assert compiler
    libfile=scratch/"tomar-six-lift-guard.so"
    subprocess.run([compiler,"-O3","-std=c++17","-ftrapv","-shared","-fPIC",
                    f"-fmodules-cache-path={cache}","-I",str(scratch),
                    str(HERE/"tomar_six_083730_guard.cpp"),"-o",str(libfile)],check=True)
    lib=ctypes.CDLL(str(libfile))
    v=lambda d:np.ctypeslib.ndpointer(dtype=d,ndim=1,flags="C_CONTIGUOUS")
    I64,I32=v(np.int64),v(np.int32)
    lib.lift_guard_run.argtypes=[I64,I64,I64,I32,I32,ctypes.c_int64,ctypes.c_int,ctypes.c_int,
          I64,I64,I64,I64,I64,ctypes.c_int,ctypes.c_int,ctypes.c_int,I64]
    tree=SourceTree(a.tree)
    T={n:load_table(getattr(a,f"table{n}"),tree,previous["tables"][str(n)]) for n in (32,64)}
    S=DirectionProof(json.loads(raw),previous["angles"])

    def scaled(values):
        x=[FINE*F(v) for v in values]
        assert all(t.denominator==1 for t in x)
        return np.array(list(map(int,x)),dtype=np.int64)

    root=scaled(tree.root["hi"])
    target=F(spec["target"])*2*Q*Q
    assert target.denominator>1

    def run(row,box,support):
        part,pos,d,mode,n,old,new,depth,tried=row
        lo,hi=(scaled(v) for v in zip(*box,strict=True))
        out=np.zeros(30,dtype=np.int64)
        result=lib.lift_guard_run(lo,hi,root,*T[d],ceil(target),depth+1,n,*support,
                                 len(support[1]),len(support[3])-1,mode,out)
        assert result==0 and out[:7].tolist()==[1,n,old,new,depth,0,tried],(part,pos,mode,out[:7])
        assert n==2*(old+new)-1
        return row

    rows=[row for group in spec["comparison"].values() for row in group] if a.compare_caps else spec["cases"]
    todo=[row for row in rows if not a.parts or row[0] in a.parts]
    assert todo
    done=[]
    with ThreadPoolExecutor(a.workers) as pool:
        queued=[]
        for row in todo:
            box=tree.leaf(*row[:2])
            queued.append(pool.submit(run,row,box,S.presolve(box)))
            if len(queued)>=a.workers*2:
                done.append(queued.pop(0).result())
        done.extend(f.result() for f in queued)
    assert done==todo
    show={"target":spec["target"],"input_sha256":digest(source),
          "all_541":len(done)==541 and not a.compare_caps,
          "covered":len(done),
          "old501_covered":sum(tuple(q[:2]) in K for q in done),
          "outer40_covered":sum(tuple(q[:2]) in B for q in done),
          "S_used_leaves":sum(q[6]>0 for q in done),"calls":sum(q[4] for q in done),
          "old_terminals":sum(q[5] for q in done),"S_terminals":sum(q[6] for q in done),
          "case_rows_sha256":digest(json.dumps(done,separators=(",",":")).encode()),
          "compiler":clang,"python":platform.python_version(),"python_flint":flint.__version__}
    if a.compare_caps:show["comparison"]=spec["comparison"]
    print(json.dumps(show,indent=2))

if __name__=="__main__":
    main()
