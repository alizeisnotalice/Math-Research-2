"""Decode authenticated preorder leaves into original integer boxes for independent witness ranking."""
import argparse
import ctypes
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
import numpy as np

SOURCE = r"""
#include <algorithm>
#include <cstdint>
using L=std::int64_t;
struct Decoder {
    const std::uint8_t* code; int n, at=0;
    const std::int32_t* pos; int total, found=0;
    L *lower,*upper;
    bool walk(L* lo,L* hi,int depth) {
        if(depth>256 || at>=n)return false;
        int index=at++;
        int axis=code[index];
        if(axis>=253) {
            if(found<total && index==pos[found]) {
                std::copy(lo,lo+11,lower+11*found);
                std::copy(hi,hi+11,upper+11*found);
                ++found;
            }
            return true;
        }
        if(axis>=11 || ((lo[axis]+hi[axis])%2))return false;
        auto a=lo[axis],b=hi[axis],m=(a+b)/2;
        hi[axis]=m;
        if(!walk(lo,hi,depth+1))return false;
        hi[axis]=b;lo[axis]=m;
        if(!walk(lo,hi,depth+1))return false;
        lo[axis]=a;
        return true;
    }
};
extern "C" int decode(const std::uint8_t* code,int n,const std::int32_t* pos,int count,L* lo,L* hi,L* lower,L* upper){
    Decoder r{code,n,0,pos,count,0,lower,upper};
    return r.walk(lo,hi,0)&&r.at==n&&r.found==count ? 0:1;
}
"""
def digest(data):return sha256(data).hexdigest()
def main():
    if not __debug__:raise SystemExit("Assertions required.")
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ("cert","seeds","seed_json","work","cache"):ap.add_argument(name,type=Path)
    args=ap.parse_args()
    sys.path.insert(0,str(Path(__file__).resolve().parent))
    from tomar_leaf_box import leaf_box
    work=args.work
    work.mkdir(parents=True,exist_ok=True)
    src=work/"decode.cpp";src.write_text(SOURCE)
    lib=work/"decode.so"
    subprocess.run(["clang++","-std=c++17","-O3","-ftrapv","-shared","-fPIC",
                    f"-fmodules-cache-path={args.cache}",src,"-o",lib],check=True)
    fn=ctypes.CDLL(str(lib.resolve())).decode
    vec=lambda dtype:np.ctypeslib.ndpointer(dtype=dtype,ndim=1,flags="C_CONTIGUOUS")
    i=ctypes.c_int;v=vec(np.int64)
    fn.argtypes=[vec(np.uint8),i,vec(np.int32),i,v,v,v,v]
    seed=json.loads(args.seed_json.read_bytes())
    raw=(args.cert/"SHA256SUMS").read_bytes()
    assert digest(raw)==seed["tree_manifest_sha256"]=="56cbe5dfdadc6b5aa07cf6ec9acf4f637a3c4fa78db54df0ae2d78916d292410"
    pins={n.removeprefix("./"):h for h,n in (s.split() for s in raw.decode().splitlines())}
    out={n:[] for n in ("part","pos","lo","hi","full_floor")}
    checks=0; scale=2**44
    for part,*row in seed["rows"]:
        stem=f"part{part:04}"
        path=args.cert/(stem+".json");codes_path=args.cert/(stem+".npy")
        assert digest(path.read_bytes())==pins[path.name] and digest(codes_path.read_bytes())==pins[codes_path.name]
        header=json.loads(path.read_bytes());assert header["part"]==part and header["split0"]==[2]*11
        lohi=[[Fraction(t)*scale for t in header[n]] for n in ("lo","hi")]
        assert all(x.denominator==1 for r in lohi for x in r)
        base_lo,base_hi=np.array(lohi,dtype=np.int64)
        assert np.all((base_lo+base_hi)%2==0)
        mid=(base_lo+base_hi)//2;bits=np.array(np.unravel_index(part,(2,)*11),bool)
        blo=np.where(bits,mid,base_lo);bhi=np.where(bits,base_hi,mid)
        with np.load(args.seeds/(stem+".npz"),allow_pickle=False) as z: ids=z["pos"];score=z["score"]
        assert ids.dtype==np.int32 and score.dtype==np.int64 and len(ids)==len(score)==row[3]
        assert digest(ids.tobytes())==row[-2] and digest(score.tobytes())==row[-1]
        codes=np.load(codes_path,allow_pickle=False)
        assert len(ids)==0 or np.all(np.diff(ids)>0) and np.all(codes[ids]>=253)
        lower=np.zeros((len(ids),11),np.int64);upper=np.empty_like(lower)
        assert fn(codes,len(codes),ids,len(ids),blo,bhi,lower.ravel(),upper.ravel())==0
        if part%17==0 and len(ids):
            prefix=np.r_[1,1+np.cumsum(np.where(codes<11,1,-1))]
            k=len(ids)//2
            c,w=leaf_box(prefix,codes,part,int(ids[k]),header)
            assert all((Fraction(x)-Fraction(y))*scale==z for x,y,z in zip(c,w,lower[k]))
            assert all((Fraction(x)+Fraction(y))*scale==z for x,y,z in zip(c,w,upper[k]))
            checks+=1
        for key,arr in (("part",np.full(len(ids),part,np.int32)),("pos",ids),("lo",lower),("hi",upper),("full_floor",score)):
            out[key].append(arr)
    arrays={k:np.concatenate(v) for k,v in out.items()}
    assert len(arrays["part"])==seed["remaining"] and np.all(arrays["lo"]<arrays["hi"])
    path=work/"candidate-endpoints.npz";np.savez_compressed(path,**arrays)
    rec=dict(target=seed["target"],coordinate_denominator=scale,area_denominator=2**57,
             dimensions=["t3x","t3y","rho5","t5x","t5y","rho7","t7x","t7y","rho9","t9x","t9y"],
             rows=seed["remaining"],original_box=True,order="part then original preorder node (both zero-based)",
             score_field="full_floor",threshold_area_integer=seed["target_integer"],
             independent_fraction_leaf_checks=checks,seed_sha256=digest(args.seed_json.read_bytes()),
             source_sha256=digest(Path(__file__).read_bytes()),file_sha256=digest(path.read_bytes()))
    (work/"candidate-endpoints.json").write_text(json.dumps(rec,indent=2)+"\n")
    print(json.dumps(rec,separators=(",",":")))
if __name__=="__main__":main()
