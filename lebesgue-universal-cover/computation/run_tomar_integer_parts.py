"""Resume atomic independent integer fan receipts for the 2,048 original 11D root cells."""
import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
import ctypes,hashlib,json,os,platform,sys,time,subprocess,tempfile
from fractions import Fraction as F
from pathlib import Path
import numpy as np

p=argparse.ArgumentParser(description=__doc__)
for key in ("certdir","outdir","library","table"):p.add_argument(key,type=Path)
p.add_argument("--target",default="0.83661")
p.add_argument("--parts",default=None)
p.add_argument("--workers",type=int,default=4)
p.add_argument("--depth",type=int,default=30)
p.add_argument("--calls",type=int,default=4095)
p.add_argument("--fresh",action="store_true",help="Recompute matching local part receipts, including their new arithmetic")
a=p.parse_args();Q=2**44;q=2**28
def digest(path):
    with path.open("rb") as f:return hashlib.file_digest(f,"sha256").hexdigest()
sys.path.insert(0,str(Path(__file__).resolve().parent))
from verify_tomar_fan_table import verified_table,PINNED_MANIFEST
first_header=json.loads((a.certdir/"part0000.json").read_bytes())
table,counts=verified_table(a.table,a.certdir,first_header)
pinned={name.removeprefix("./"):sha for sha,name in
        (line.split() for line in (a.certdir/"SHA256SUMS").read_text().splitlines())}
assert digest(a.certdir/"part0000.json")==pinned["part0000.json"]
# The supplied path names an OUTPUT binary. Rebuild from this source on every run.
# An unrelated existing library must never be accepted by recording its digest.
cpp=Path(__file__).with_name("tomar_integer_fan.cpp").resolve()
compiler=Path("/usr/bin/clang++")
cache=Path("/Volumes/data/CodexTemp/cache/clang")
assert cache.parent.parent.is_dir()
cache.mkdir(parents=True,exist_ok=True)
a.library.parent.mkdir(parents=True,exist_ok=True)
flags=["-O3","-ftrapv","-std=c++17","-fmodules-cache-path="+str(cache),
       "-shared","-fPIC","-Wl,-install_name,@rpath/libtomar_fan.dylib"]
with tempfile.TemporaryDirectory(dir=a.library.parent) as tmp:
    binary=Path(tmp)/"libtomar_fan.dylib"
    subprocess.run([str(compiler),*flags,str(cpp),"-o",str(binary)],check=True)
    compiled_sha=digest(binary)
    if a.library.is_file():
        if digest(a.library)!=compiled_sha:
            raise ValueError("the specified library differs from the reproducible local build")
    else:
        binary.replace(a.library)
compiler_version=subprocess.check_output([str(compiler),"--version"],text=True).splitlines()[0]
lib=ctypes.CDLL(str(a.library.resolve()))
vec=lambda t: np.ctypeslib.ndpointer(dtype=t,ndim=1,flags="C_CONTIGUOUS")
lib.adaptive_trim.argtypes=[vec(np.uint8),ctypes.c_int,vec(np.int64),vec(np.int64),vec(np.int64),
    vec(np.uint8),vec(np.uint8),vec(np.int64),vec(np.int32),vec(np.int32),ctypes.c_int64,ctypes.c_int,ctypes.c_int]
lib.adaptive_trim.restype=ctypes.c_int
threshold=2*F(a.target)*q*q;goal=-(-threshold.numerator//threshold.denominator)
commitment={"target":a.target,"max_depth":a.depth,"max_calls":a.calls,"scale":q,"table_grid":1024,
    "table_sha256":digest(a.table),"library_sha256":compiled_sha,"compiler_sha256":digest(compiler),"compiler_version":compiler_version,"compiler_flags":flags,
    "input_manifest_sha256":PINNED_MANIFEST,"generator_sha256":digest(Path(__file__).with_name("tomar_fan_table.py")),
    "verifier_sha256":digest(Path(__file__).with_name("verify_tomar_fan_table.py")),"driver_sha256":digest(Path(__file__)),
    "cpp_source_sha256":digest(Path(__file__).with_name("tomar_integer_fan.cpp"))}
a.outdir.mkdir(parents=True,exist_ok=True)
keys=range(2048) if a.parts is None else sorted(set(map(int,a.parts.split(","))))
def process(part):
    assert 0<=part<2048
    stem=f"part{part:04}";meta_path=a.certdir/(stem+".json");codes_path=a.certdir/(stem+".npy")
    basis={**commitment,"part":part,"header_sha256":digest(meta_path),"codes_sha256":digest(codes_path),
           "prior_unresolved_sha256":None,"initial_mask":"all original nodes with code 253, 254 or 255"}
    assert basis["header_sha256"]==pinned[meta_path.name] and basis["codes_sha256"]==pinned[codes_path.name]
    dest=a.outdir/(stem+".json");out=a.outdir/(stem+".npz")
    if dest.exists() and out.is_file():
        r=json.loads(dest.read_bytes())
        if not (all(r.get(k)==v for k,v in basis.items()) and digest(out)==r["remaining_sha256"]):
            raise ValueError("existing receipt differs at "+stem)
        cached_codes=np.load(codes_path,allow_pickle=False)
        cached_pos=np.load(out,allow_pickle=False)["pos"]
        original_count=int(np.count_nonzero(cached_codes>=253))
        assert original_count==r["all_original"]==r["tested"]==r["closed"]+r["remaining"]
        assert cached_pos.dtype==np.int32 and np.array_equal(cached_pos,np.unique(cached_pos))
        assert ((0<=cached_pos)&(cached_pos<len(cached_codes))).all() and (cached_codes[cached_pos]>=253).all()
        assert len(cached_pos)==r["remaining"]
        if not a.fresh:return r
    begin=time.monotonic();meta=json.loads(meta_path.read_bytes());codes=np.load(codes_path)
    assert meta["ns"]==[3,5,7,9] and meta["split0"]==[2]*11 and meta["nparts"]==2048 and meta["part"]==part
    assert meta["lo"]==first_header["lo"] and meta["hi"]==first_header["hi"]
    assert codes.dtype==np.uint8 and codes.ndim==1
    old=(codes>=253).astype(np.uint8)
    num_all=int(old.sum())
    lo=[];hi=[];original=[]
    for x,y,bit in zip(meta["lo"],meta["hi"],np.unravel_index(part,(2,)*11)):
        l,h=F(x)*Q,F(y)*Q;assert l.denominator==h.denominator==1 and -Q<l<h<Q
        original.append(int(h));l,h=((l+h)/2,h) if bit else (l,(l+h)/2)
        assert l.denominator==h.denominator==1;lo.append(int(l));hi.append(int(h))
    mask=np.zeros_like(old);s=np.zeros(6,dtype=np.int64)
    rc=lib.adaptive_trim(codes,len(codes),np.array(lo),np.array(hi),np.array(original),old,mask,s,
        table.ravel(),counts,goal,a.depth,a.calls)
    assert rc==0,(part,rc)
    pos=np.flatnonzero(mask).astype(np.int32)
    assert np.all(mask<=old) and s[2]+s[3]==old.sum() and s[3]==len(pos)
    tmp=out.with_suffix(".tmp.npz");np.savez_compressed(tmp,pos=pos);tmp.replace(out)
    record={**basis,"all_original":num_all,"tested":int(old.sum()),"closed":int(s[2]),"remaining":int(s[3]),
        "checked_boxes":int(s[0]),"maximum_depth_used":int(s[1]),"closed_subboxes_from_successes":int(s[4]),
        "budget_limited":int(s[5]),"remaining_sha256":digest(out),
        "remaining_codes":{str(k):int(np.count_nonzero(codes[pos]==k)) for k in [253,254,255]},
        "wall_seconds":round(time.monotonic()-begin,3)}
    tmp=dest.with_suffix(".partial");tmp.write_text(json.dumps(record,indent=2)+"\n");tmp.replace(dest)
    print(f"part {part:04}, {record['closed']}/{record['tested']} closed",file=sys.stderr,flush=True)
    return record
keys=sorted(keys,key=lambda k:(a.certdir/f"part{k:04}.npy").stat().st_size,reverse=True)
with ThreadPoolExecutor(max_workers=a.workers) as pool:
    futures=[pool.submit(process,k) for k in keys]
    rows=sorted((f.result() for f in as_completed(futures)),key=lambda r:r["part"])
print(json.dumps({**commitment,"platform":platform.platform(),"python":platform.python_version(),
    "parts":len(rows),"tested":sum(r["tested"] for r in rows),"closed":sum(r["closed"] for r in rows),
    "remaining":sum(r["remaining"] for r in rows),"boxes":sum(r["checked_boxes"] for r in rows),
    "rows":[[r["part"],r["tested"],r["closed"],r["remaining"],r["checked_boxes"],r["wall_seconds"]] for r in rows]},separators=(",",":")))
