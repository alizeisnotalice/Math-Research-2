"""Verify the 30-point extension on all residual original leaves and on two stricter blocks."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import ctypes
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
import platform
import subprocess
import sys

import flint
from flint import arb,ctx
import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from certify_finite_point_leaf import digest, ar, upper_q, point_tables, FINE,Q,U,TABLE_SHA
from verify_tomar_fan_table import verified_table,PINNED_MANIFEST

DENSE_SHA="0913f79e54bfcb804c8c4a4c0bbd7962d42d461591175ab037739b3562c6016d"


def run():
    parser=argparse.ArgumentParser(description=__doc__)
    for k in ("cert","table32","table64","work","cache"):
        parser.add_argument(k,type=Path)
    parser.add_argument("--workers",type=int,default=1)
    parser.add_argument("--parts",help="verify only the stated roots, for example 180")
    parser.add_argument("--extra-only",action="store_true",help="verify the two additional .83730 old-low leaves")
    parser.add_argument("--replay-five",action="store_true",help="also recompute the outside-of-501 five-body proof")
    parser.add_argument("--lift",type=Path,help="higher-threshold terminal witness")
    parser.add_argument("--outer",type=Path,help="higher-threshold five-body outer leaf index")
    a=parser.parse_args()
    assert __debug__ and 1<=a.workers<=32
    work=a.work.resolve()
    work.mkdir(exist_ok=False)
    cache=a.cache.resolve()
    cache.mkdir(parents=True,exist_ok=True)
    ctx.prec=160
    inp=HERE/"finite-point-global-input.json"
    rec=json.loads(inp.read_bytes())
    old_idx=HERE/"five-global-localization.json"
    old_receipt=HERE/"five-global-localization-replay.json"
    cap=F(rec["target"])
    old=json.loads(old_idx.read_bytes())
    five=json.loads(old_receipt.read_bytes())
    assert five["complete"] and five["expected_sha256"]==digest(old_idx) and five["parts"]==2048
    assert five["totals"]==old["totals"] and F(old["target"])==cap==F(".837295006222")
    assert a.cert.is_dir() and digest(a.cert/"SHA256SUMS")==PINNED_MANIFEST
    manifest={n.removeprefix("./"):v for v,n in
              (r.split() for r in (a.cert/"SHA256SUMS").read_text().splitlines())}
    if a.replay_five:
        assert not a.extra_only and a.parts is None
        subprocess.run([sys.executable,"-I","-B",HERE/"replay_five_global_localization.py",
                        a.cert,work/"five-exterior",cache,"--workers",str(a.workers),"--table",a.table32],check=True)
    keys=sorted((int(k),idx) for k,r in old["retained"].items() for idx in r)
    assert len(keys)==501
    rows=rec["rows"]
    assert [(r[0],r[1]) for r in rows]==keys
    assert (a.lift is None)==(a.outer is None)
    if a.lift:
        assert not a.extra_only
        lifted=json.loads(a.lift.read_bytes())
        outer=json.loads(a.outer.read_bytes())
        assert digest(inp)==lifted["prior_bank_sha256"] and digest(a.outer)==lifted["outer_sha256"]
        assert outer["tree_manifest_sha256"]==PINNED_MANIFEST
        cap=F(lifted["target"])
        assert cap>F(rec["target"]) and F(outer["target"])==cap
        extra_keys={(int(p),i) for p,x in outer["retained"].items() for i in x}
        assert not extra_keys.intersection(keys)
        keys=sorted(set(keys)|extra_keys)
        rows=lifted["rows"]
        assert [(r[0],r[1]) for r in rows]==keys

    headfile=a.cert/"part0000.json"
    assert digest(headfile)==manifest[headfile.name]
    head=json.loads(headfile.read_bytes())
    root_lo,root_hi=(np.array([int(F(x)*FINE) for x in head[k]],np.int64) for k in ("lo","hi"))
    assert all(F(x)*FINE==y for k,vals in (("lo",root_lo),("hi",root_hi)) for x,y in zip(head[k],vals))
    tables={}
    for k,path,want in ((32,a.table32,TABLE_SHA),(64,a.table64,DENSE_SHA)):
        assert digest(path)==want
        pts,cnt=verified_table(path,a.cert,head)
        tables[k]=pts.ravel(),cnt

    shape_path=HERE/"finite-point-leaf-input.json"
    shape=json.loads(shape_path.read_bytes())
    assert len(shape["points"])==30
    shape["cells"]=rec["cells"]
    dia2,ptr,dual=point_tables(shape)
    nums=shape["directions"]
    vecs=[(ar(F(x,U)),ar(F(y,U)),ar(F(x*x+y*y,U*U)).sqrt()) for x,y in nums]
    disk=[upper_q(length/2) for _,_,length in vecs]
    src=HERE/"finite_point_leaf_guard.cpp"
    binary=work/"point-guard.dylib"
    flags=["-O3","-std=c++17","-shared","-fPIC","-fsanitize=signed-integer-overflow",
           "-fsanitize-undefined-trap-on-error","-fmodules-cache-path="+str(cache)]
    subprocess.run(["/usr/bin/clang++",*flags,src,"-o",binary],check=True)
    fn=ctypes.CDLL(str(binary)).dynamic_guard_run
    v=lambda typ:np.ctypeslib.ndpointer(dtype=typ,ndim=1,flags="C_CONTIGUOUS")
    i64,i32=v(np.int64),v(np.int32)
    fn.argtypes=[i64,i64,i64,i32,i32,ctypes.c_int64,ctypes.c_int,ctypes.c_int,
                 i64,i64,i64,i64,i64,ctypes.c_int,ctypes.c_int,i64]
    fn.restype=ctypes.c_int
    ints=lambda seq:np.ascontiguousarray(seq,dtype=np.int64).ravel()
    ns,ds,ptr,dual=map(ints,(nums,disk,ptr,dual))

    @lru_cache(None)
    def root(p):
        header=a.cert/f"part{p:04}.json"
        file=a.cert/f"part{p:04}.npy"
        assert digest(header)==manifest[header.name] and digest(file)==manifest[file.name]
        h=json.loads(header.read_bytes())
        assert h["part"]==p and h["nparts"]==2048 and h["split0"]==[2]*11 and h["ns"]==[3,5,7,9]
        assert h["lo"]==head["lo"] and h["hi"]==head["hi"]
        codes=np.load(file,allow_pickle=False)
        assert codes.dtype==np.uint8 and codes.ndim==1 and len(codes)==h["boxes"]
        stop=[0]*len(codes)
        for i,k in reversed(list(enumerate(codes.tolist()))):
            if k>=253:stop[i]=i+1
            else:
                assert 0<=k<11 and i+1<len(codes) and stop[i+1]<len(codes)
                stop[i]=stop[stop[i+1]]
        assert stop[0]==len(codes)
        return codes,stop

    def leaf(p,pos):
        codes,stop=root(p)
        lo,hi=root_lo.copy(),root_hi.copy()
        mid=(lo+hi)//2
        bits=np.array(np.unravel_index(p,(2,)*11),bool)
        lo[bits],hi[~bits]=mid[bits],mid[~bits]
        node=0
        while node!=pos:
            dim=int(codes[node])
            assert dim<11 and node<pos<stop[node] and (lo[dim]+hi[dim])%2==0
            mid=(lo[dim]+hi[dim])//2
            if pos<stop[node+1]:
                hi[dim]=mid
                node+=1
            else:
                lo[dim]=mid
                node=stop[node+1]
        assert codes[node]>=253
        return lo,hi

    @lru_cache(None)
    def old_nominal(n,turn):
        halfcos=(arb.pi()/(2*n)).cos()
        radius=1/(2*halfcos)
        assert radius<ar(F(3,5))
        vs=[(ar(F(turn,FINE))+2*arb.pi()*i/n).sin_cos() for i in range(n)]
        result=[]
        for x,y,length in vecs:
            u=[]
            for s,c in vs:
                d=c*x+s*y
                u.append(radius*d)
                if not d>-length*halfcos:u.append(length+radius*d)
            result.append(max(upper_q(z) for z in u))
        return result

    def execute(row,target):
        p,pos,k,budget,depth,vis,old_count,new_count,maxdepth=row
        lo,hi=leaf(p,pos)
        ends=[]
        for n,axis in ((3,None),(5,2),(7,5),(9,8)):
            low,high=(0,0) if axis is None else (int(lo[axis]),int(hi[axis]))
            err=Q*F(3,5)*F(high-low,FINE)**2*F(U+2,U)/8
            u=-(-err.numerator//err.denominator)
            for ang in (low,high):
                ends.extend(z+u for z in old_nominal(n,ang))
        ret=np.zeros(29,np.int64)
        pts,cnt=tables[k]
        threshold=F(target)*2*Q*Q
        goal=threshold.numerator//threshold.denominator+1
        assert fn(lo,hi,root_hi,pts,cnt,goal,depth,budget,ns,ds,ints(ends),ptr,dual,
                  len(ds),len(ptr)-1,ret)==0
        assert ret[:6].tolist()==[1,vis,old_count,new_count,maxdepth,0],(p,pos,ret[:7].tolist(),row)
        assert vis==2*(old_count+new_count)-1
        return ret[1:4].tolist()

    wanted=set(int(x) for x in a.parts.split(",")) if a.parts else set(p for p,_ in keys)
    assert wanted<=set(p for p,_ in keys)
    scores=[]
    if not a.extra_only:
        use=[r for r in rows if r[0] in wanted]
        with ThreadPoolExecutor(a.workers) as pool:
            scores=list(pool.map(lambda r:execute(r,cap),use))
    from certify_five_body_contact_stationary import model,area,verify_topology
    extra=[]
    for r in ([] if a.lift else rec["extra_blocks"]):
        if r["part"] not in wanted:continue
        ctx.prec=192
        alpha,beta=ar(F(r["alpha"])),ar(F(r["beta"]))
        topology=verify_topology(alpha,beta)
        pose=model(alpha,beta)[1]
        lo,hi=leaf(r["part"],r["pos"])
        assert all(ar(F(int(l),FINE))<x and x<ar(F(int(h),FINE)) for l,h,x in zip(lo,hi,pose))
        value=area(alpha,beta)
        assert value<ar(F(".83730")) and r["pos"]!=23955
        ctx.prec=160
        row=[r["part"],r["pos"],32,4000000,99,r["calls"],r["old"],r["new"],r["depth"]]
        execute(row,F(".83730"))
        extra.append(dict(part=r["part"],pos=r["pos"],old_area=str(value),topology=topology,visits=r["calls"]))
    out=dict(complete=not a.extra_only and a.parts is None,old_root_replayed=a.replay_five,
             target=str(cap),original_leaves=len(scores),successful_tree_totals=np.sum(scores,axis=0,dtype=np.int64).tolist() if scores else [],
             nontrivial_point_leaves=sum(v[2]>0 for v in scores),extra_blocks=extra,points=30,diameter_squared_numerator=dia2,
             directions=len(ds),cells=len(ptr)-1,balanced_triples=len(dual)//8,prior_localization=digest(old_idx),
             geometry_sha256=digest(shape_path),bank_and_rows=digest(inp),kernel_sha256=digest(src),binary_sha256=digest(binary),
             driver_sha256=digest(Path(__file__)),table32=TABLE_SHA,table64=DENSE_SHA,tree_manifest=PINNED_MANIFEST,
             compiler_flags=flags,python=platform.python_version(),flint=flint.__version__)
    if a.lift:
        out.update(lift_sha256=digest(a.lift),outer_sha256=digest(a.outer),
                   prior_score_complement="see the independently audited higher-threshold five-body exterior")
    (work/"replay.json").write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps(out,indent=2))


if __name__=="__main__":
    run()
