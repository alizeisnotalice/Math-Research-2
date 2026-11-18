"""Apply fresh five-body prefilters at .83730 to the independently selected outer leaves."""
import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
import ctypes
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import time
import numpy as np

TARGET=F(".83730")
FINE=2**44
Q=2**28
S=2*Q*Q
PINNED_MANIFEST="56cbe5dfdadc6b5aa07cf6ec9acf4f637a3c4fa78db54df0ae2d78916d292410"
TABLE_HASH="fe5b816060d713642dc44e3f8981f4ae5d1b2b1c1e271b7dff75f031dd75d87c"
STAGES=((2,127),(1,1023),(1,8191),(1,65535))

def sha(x):
    return sha256(x.read_bytes()).hexdigest()

def ceil(q):
    return -(-q.numerator//q.denominator)

def run():
    p=argparse.ArgumentParser(description=__doc__)
    for x in ("tree","outer","old501","table","radial","kernel","work"):p.add_argument(x,type=Path)
    p.add_argument("--workers",type=int,default=3)
    p.add_argument("--score-audit",type=Path,required=True)
    p.add_argument("--expected",type=Path)
    p.add_argument("--clang-cache",type=Path)
    p.add_argument("--parts",nargs="+",type=int)
    a=p.parse_args()
    a.work.mkdir(parents=True,exist_ok=False)
    C=a.tree
    verified=json.loads(a.score_audit.read_bytes())
    assert F(verified["target"])==TARGET
    assert verified["tree_manifest_sha256"]==PINNED_MANIFEST
    assert verified["totals"]["remaining_outer"]==1026625
    audited={r[0]:r for r in verified["rows"]}
    assert len(audited)==1024
    assert sha(C/"SHA256SUMS")==PINNED_MANIFEST
    M={name.removeprefix("./"):h for h,name in (s.split() for s in (C/"SHA256SUMS").read_text().splitlines())}
    assert sha(C/"part0000.json")==M["part0000.json"]
    base=json.loads((C/"part0000.json").read_bytes())
    lo0=np.array([int(F(q)*FINE) for q in base["lo"]],dtype=np.int64)
    hi0=np.array([int(F(q)*FINE) for q in base["hi"]],dtype=np.int64)
    mid0=(lo0+hi0)//2
    assert np.all(2*mid0==lo0+hi0)
    radial=json.loads(a.radial.read_bytes())["rows"]
    limits=[F(r["radius"]) for r in radial]
    assert [r["n"] for r in radial]==[3,5,7,9]
    assert all(F(x["uniform_twice_area_numerator"],S)>TARGET for x in radial)
    cpp=a.kernel.read_text()
    assert "r[4]={"+",".join(str(ceil(x*FINE)) for x in limits)+"};" in cpp
    matrix={}
    assert sha(a.table)==TABLE_HASH
    with np.load(a.table,allow_pickle=False) as data:
        table,counts=data["points"],data["counts"]
    for div,_ in STAGES:
        offs=np.r_[0,np.cumsum(counts)]
        jj=[j for l,r in zip(offs[:-1],offs[1:],strict=True) for j in range(l,r,div)]
        matrix[div]=np.ascontiguousarray(table[:,jj,:]).ravel(),np.array([(int(q)+div-1)//div for q in counts],np.int32)
    binary=a.work/"outer-five.dylib"
    subprocess.run(["/usr/bin/clang++","-std=c++17","-O3","-ftrapv","-shared","-fPIC",
        f"-fmodules-cache-path={a.clang_cache or a.work/'clang-cache'}",str(a.kernel),"-o",str(binary)],check=True)
    lib=ctypes.CDLL(str(binary.resolve()))
    v=lambda t:np.ctypeslib.ndpointer(dtype=t,ndim=1,flags="C_CONTIGUOUS")
    U,I,J=v(np.uint8),v(np.int64),v(np.int32)
    lib.adaptive_trim.argtypes=[U,ctypes.c_int,I,I,I,U,U,I,J,J,ctypes.c_int64,ctypes.c_int,ctypes.c_int]
    previous=json.loads(a.old501.read_bytes())
    if "retained" in previous:
        previous=previous["retained"]
    assert sum(map(len,previous.values()))==501
    files=list(a.outer.glob("part*.npz"))
    assert len(files)==1024
    if a.parts:
        files=[f for f in files if int(f.stem[4:]) in a.parts]
        assert len(files)==len(set(a.parts))
    goal=ceil(TARGET*S)
    def one(src):
        name=src.stem
        part=int(name[4:])
        path=C/(name+".npy")
        hdr=C/(name+".json")
        assert sha(path)==M[path.name] and sha(hdr)==M[hdr.name]
        meta=json.loads(hdr.read_bytes())
        assert meta["part"]==part and all(meta[k]==base[k] for k in ("ns","lo","hi","split0","nparts"))
        codes=np.load(path,allow_pickle=False)
        with np.load(src,allow_pickle=False) as q:
            ids,score,inside=q["pos"],q["score"],q["inside501"]
            assert ids.ndim==score.ndim==inside.ndim==1 and len(ids)==len(score)==len(inside)
            assert np.array_equal(ids,np.unique(ids)) and ((score<=2000)).all()
            row=audited[part]
            assert sha256(ids.tobytes()).hexdigest()==row[-2]
            assert sha256(score.tobytes()).hexdigest()==row[-1]
            assert ids[inside].tolist()==previous.get(str(part),[])
            other=ids[~inside]
        assert len(other)==audited[part][6]
        assert np.all((0<=other)&(other<len(codes))) and (codes[other]>=253).all()
        mask=np.zeros(len(codes),dtype=np.uint8)
        mask[other]=1
        bits=np.array(np.unravel_index(part,(2,)*11),dtype=bool)
        lo,hi=lo0.copy(),hi0.copy()
        lo[bits],hi[~bits]=mid0[bits],mid0[~bits]
        seq=[len(other)]
        vis=[]
        for div,cap in STAGES:
            next_=np.zeros_like(mask)
            stats=np.zeros(6,dtype=np.int64)
            if mask.any():
                assert lib.adaptive_trim(codes,len(codes),lo,hi,hi0,mask,next_,stats,*matrix[div],goal,48,cap)==0
            assert (next_<=mask).all() and next_.sum()+stats[2]==mask.sum()
            mask=next_
            seq.append(int(mask.sum()))
            vis.append(int(stats[0]))
        final=np.flatnonzero(mask).tolist()
        return {"part":part,"input_sha256":sha(src),"raw_codes_sha256":M[path.name],
                "raw_header_sha256":M[hdr.name],"counts":seq,"visits":vis,"retained":final}

    rows=[]
    begin=time.monotonic()
    with ThreadPoolExecutor(max_workers=a.workers) as pool,(a.work/"parts.jsonl").open("w") as out:
        for f in as_completed([pool.submit(one,s) for s in sorted(files)]):
            x=f.result()
            rows.append(x);out.write(json.dumps(x,separators=(",",":"))+"\n");out.flush()
            if len(rows)%64==0:
                print(len(rows),[sum(r["counts"][i] for r in rows) for i in range(5)],
                      round(time.monotonic()-begin,2),flush=True)
    rows.sort(key=lambda v:v["part"])
    totals=[sum(r["counts"][i] for r in rows) for i in range(5)]
    assert totals[0]==sum(audited[r["part"]][6] for r in rows)
    output={"complete":a.parts is None,"score_audit_sha256":sha(a.score_audit),"target":str(TARGET),"tree_manifest_sha256":PINNED_MANIFEST,"kernel_sha256":sha(a.kernel),
            "validated_old_table_sha256":TABLE_HASH,"radial_archive_sha256":sha(a.radial),"old501_sha256":sha(a.old501),
            "source_outside501":totals[0],"stages":STAGES,"totals":totals,
            "stage_node_visits":[sum(r["visits"][i] for r in rows) for i in range(4)],
            "parent_receipts_sha256":sha(a.work/"parts.jsonl"),
            "retained":{str(r["part"]):r["retained"] for r in rows if r["retained"]}}
    if a.expected:
        expected=json.loads(a.expected.read_bytes())
        assert F(expected["target"])==TARGET
        assert output["retained"]=={
            key:value for key,value in expected["retained"].items() if not a.parts or int(key) in a.parts}
        if not a.parts:
            assert output["totals"]==expected["totals"] and output["stage_node_visits"]==expected["stage_node_visits"]
    (a.work/"outer-retained.json").write_text(json.dumps(output,indent=2)+"\n")
    print("done",totals,flush=True)
if __name__=="__main__":run()
