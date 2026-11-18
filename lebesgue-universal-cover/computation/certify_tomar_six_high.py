"""Replay the 501 whole-leaf certificates at the high five-body localization threshold."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import ctypes
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
import platform
import shutil
import subprocess

from flint import arb, ctx
import flint
import numpy as np

from audit_five_harmonic_tree import read_codes, split
from certify_large_support import curvature, normal
from certify_nine import ar, PI_UPPER
from verify_tomar_fan_table import verified_table, PINNED_MANIFEST

HERE = Path(__file__).resolve().parent
U, Q, FINE = 2**22, 2**28, 2**44


def digest(data):
    return sha256(data).hexdigest()


def ceil(x):
    return -(-x.numerator//x.denominator)


class SourceTree:
    def __init__(self,path):
        self.path=path
        manifest=(path/"SHA256SUMS").read_bytes()
        assert digest(manifest)==PINNED_MANIFEST
        self.pins={name.removeprefix("./"):h for h,name in (line.split() for line in manifest.decode().splitlines())}
        self.root=json.loads(self.read("part0000.json"))
        assert self.root["ns"]==[3,5,7,9] and self.root["split0"]==[2]*11

    def read(self,name):
        data=(self.path/name).read_bytes()
        assert digest(data)==self.pins[name]
        return data

    @lru_cache
    def load(self,part):
        s=f"part{part:04}"
        d=json.loads(self.read(s+".json"))
        codes,prefix=read_codes(self.read(s+".npy"))
        assert d["part"]==part and d["nparts"]==2048 and len(codes)==d["boxes"]
        assert all(d[k]==self.root[k] for k in ("lo","hi","split0","ns"))
        box=list(zip(map(F,d["lo"]),map(F,d["hi"]),strict=True))
        for k in range(11):
            box=split(box,k)[(part>>(10-k))&1]
        return box,codes,prefix

    def leaf(self,part,pos):
        box,codes,prefix=self.load(part)
        assert 0<=pos<len(codes) and codes[pos]>=253
        k=0
        while k!=pos:
            axis=codes[k]
            assert axis<11
            r=prefix.index(prefix[k],k+2)
            side=int(pos>=r)
            box=split(box,axis)[side]
            k=r if side else k+1
        return box


def outward(z):
    q=round(float(z)*Q)+2
    assert ar(F(q,Q))>z
    return q


def old_support(n,angle,vector,pi):
    c=(pi/(2*n)).cos()
    r=1/(2*c)
    x,y,length=vector
    a=[]
    for j in range(n):
        t=ar(angle)+2*pi*j/n
        d=x*t.cos()+y*t.sin()
        a.append(r*d)
        if not d> -c*length:
            a.append(length+r*d)
    return max(map(outward,a))


def new_support(vector,angle,terms):
    x,y,length=vector
    c=(x*angle.cos()+y*angle.sin())/length
    cosines=[arb(1),c]
    for _ in range(2,terms[-1][0]+1):
        cosines.append(2*c*cosines[-1]-cosines[-2])
    return length*(ar(F(1,2))+sum(ar(a)*cosines[n] for n,a in terms))


class DirectionProof:
    def __init__(self,design,choices):
        ctx.prec=160
        self.pi=arb.pi()
        assert self.pi<ar(PI_UPPER)
        terms=[(k,F(a)) for k,a in design["terms"]]
        assert curvature(terms,self.pi)==F(1,100)
        raw=[tuple(round(q*U) for q in normal(s)) for s in design["normal_seeds"]]
        self.v=[]
        for x,y in raw:
            z=x*x+y*y
            assert (U-2)**2<z<(U+2)**2
            self.v.append((ar(F(x,U)),ar(F(y,U)),ar(F(z,U*U)).sqrt()))
        self.directions=np.array(raw,dtype=np.int64).ravel()
        self.disk=np.array([outward(length/2) for _,_,length in self.v],dtype=np.int64)
        count=design["angle_cells"]
        assert len(choices)==count
        E=sum(k*k*abs(a) for k,a in terms)*F(U+2,U)*(2*PI_UPPER/count)**2/8
        points=[[new_support(v,2*self.pi*j/count,terms) for v in self.v] for j in range(count+1)]
        jumps=[0]
        rows=[]
        def det(x,y):
            return x[0]*y[1]-x[1]*y[0]
        for j,group in enumerate(choices):
            assert group
            for inds in group:
                assert len(inds)==len(set(inds))==3
                assert all(type(i) is int and 0<=i<len(raw) for i in inds)
                a,b,c=[raw[k] for k in inds]
                w=[det(b,c),det(c,a),det(a,b)]
                if max(w)<0:
                    w=[-q for q in w]
                assert min(w)>0
                assert all(sum(q*u[k] for q,u in zip(w,(a,b,c),strict=True))==0 for k in (0,1))
                t=sum(w)
                lo=[]
                for h in (j,j+1):
                    value=sum(v*points[h][k] for v,k in zip(w,inds,strict=True))/t-ar(E)
                    q=round(float(value)*Q)-2
                    assert ar(F(q,Q))<value
                    lo.append(q)
                rows.extend([*inds,*w,t,min(lo)])
            jumps.append(len(rows)//8)
        self.pointers=np.array(jumps,dtype=np.int64)
        self.triples=np.array(rows,dtype=np.int64)

    @lru_cache(maxsize=None)
    def one(self,n, angle):
        assert 1/(2*(self.pi/(2*n)).cos())<ar(F(3,5))
        return np.array([old_support(n,angle,v,self.pi) for v in self.v],dtype=np.int64)

    def presolve(self,box):
        ends=[]
        for n,index in ((3,None),(5,2),(7,5),(9,8)):
            a,b=(F(0),F(0)) if index is None else box[index]
            pad=ceil(Q*F(3,5)*(b-a)**2*F(U+2,U)/8)
            ends.extend(self.one(n,t)+pad for t in (a,b))
        return (self.directions,self.disk,np.array(ends,dtype=np.int64).ravel(),self.pointers,self.triples)


def build_kernel(scratch,spec,cache):
    source=(HERE/"tomar_grouped_scored_fan.cpp").read_bytes()
    assert digest(source)==spec["old_base_source_sha256"]
    s=source.decode().split("\nstruct Runner",1)[0]
    changes=[
        ("\n\nW up","\nW up"),
        ("    L latest_floor=0,minimum_floor=std::numeric_limits<L>::max();\n",""),
        ("L goal;int total","W last=0; L goal;int total"),
        ("L hx[5]{}","L rot_lo[5]{},rot_hi[5]{},width[5]{},hx[5]{}"),
        ("pts.clear();hull.clear();conv.clear();int off=0,ends[5]{};",
         "last=0;pts.clear();hull.clear();conv.clear();int off=0,ends[5]{};"),
        ("theta[b]=small(up(std::max(std::abs(1024*lo[k]-anchor),std::abs(1024*hi[k]-anchor)),1024*ratio));",
         """L den=1024*ratio;
                    L lower=1024*lo[k]-anchor,upper=1024*hi[k]-anchor;
                    rot_lo[b]=lower>=0?lower/den:-small(up(-lower,den));
                    rot_hi[b]=upper<=0?-((-upper)/den):small(up(upper,den));
                    width[b]=small(up(hi[k]-lo[k],2*ratio));
                    theta[b]=std::max(std::abs(rot_lo[b]),std::abs(rot_hi[b]));"""),
        ("error+=up(ab(grot[k])*theta[k],Q);",
         """// Actual interval displacement from its witness row can be one-sided.
            // |sin(a)-a| <= |a|^3/6 for every real a.
            W scalar=grot[k]*(grot[k]>=0?rot_lo[k]:rot_hi[k]);
            twice+=scalar>=0?scalar/Q:-up(-scalar,Q);
            W first=up(ab(grot[k])*theta[k],Q);
            error+=std::min(first,up(up(first*theta[k],Q)*theta[k],6*Q));"""),
        ("if(score>best){split_axis=index;best=score;}",
         "if(hi[index]>lo[index] && (hi[index]+lo[index])%2==0 && score>best){split_axis=index;best=score;}"),
        ("if(b>1)consider(k-1,W(theta[b])*(ab(grot[b])/Q+Q/100));",
         "if(b>1)consider(k-1,W(width[b])*(ab(grot[b])/Q+Q/100));"),
        ("""W bound=twice-error;
        latest_floor=bound>0?static_cast<L>(std::min<W>(bound,std::numeric_limits<L>::max())):0;
        return bound>=goal ? 1 : 0;""","last=twice-error;\n        return last>=goal ? 1 : 0;"),
    ]
    for a,b in changes:
        assert s.count(a)==1,a
        s=s.replace(a,b)
    s=s.rstrip()+"\n"
    assert digest(s.encode())==spec["exposed_header_sha256"],digest(s.encode())
    scratch.mkdir(parents=True,exist_ok=True)
    (scratch/"fan-exposed.hpp").write_text(s)
    path=scratch/"tomar-six-high-guard.so"
    clang=shutil.which("clang++")
    assert clang
    flags=["-O3","-std=c++17","-ftrapv",f"-fmodules-cache-path={cache}","-shared","-fPIC"]
    subprocess.run([clang,*flags,"-I",str(scratch),str(HERE/"tomar_six_high_guard.cpp"),"-o",str(path)],check=True)
    return ctypes.CDLL(str(path)),subprocess.check_output([clang,"--version"],text=True).splitlines()[0]


def load_table(table,tree,spec):
    points,counts=verified_table(table,tree.path,tree.root)
    stamp=sha256()
    for a in (counts,points):
        stamp.update(np.asarray(a,dtype="<i4").tobytes())
    assert stamp.hexdigest()==spec["canonical_arrays_sha256"]
    return points.ravel(),counts


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("tree","localization","table32","table64","work"):p.add_argument(k,type=Path)
    p.add_argument("--clang-cache",type=Path)
    p.add_argument("--workers",type=int,default=3)
    p.add_argument("--parts",nargs="+",type=int)
    a=p.parse_args()
    source=(HERE/"tomar-six-high-input.json").read_bytes()
    spec=json.loads(source)
    loc=json.loads(a.localization.read_bytes())
    selected=loc["retained"]
    assert F(loc["target"])==F(spec["target"]) and loc["denominator"]==FINE
    assert digest(json.dumps(selected,sort_keys=True,separators=(",",":")).encode())==spec["localization_retained_sha256"]
    cases=spec["leaves"]
    assert len(cases)==501 and [(p,i) for p,i,*_ in cases]==sorted((int(p),i) for p,L in selected.items() for i in L)
    d=(HERE/"large-support-harmonic-input.json").read_bytes()
    assert digest(d)==spec["design_sha256"]
    work=a.work.resolve()
    guard,clang=build_kernel(work,spec,a.clang_cache or work/"clang-cache")
    tree=SourceTree(a.tree)
    T={n:load_table(getattr(a,f"table{n}"),tree,spec["tables"][str(n)]) for n in (32,64)}
    S=DirectionProof(json.loads(d),spec["angles"])
    v=lambda d:np.ctypeslib.ndpointer(dtype=d,ndim=1,flags="C_CONTIGUOUS")
    I64,I32=v(np.int64),v(np.int32)
    guard.dynamic_guard_run.argtypes=[I64,I64,I64,I32,I32,ctypes.c_int64,ctypes.c_int,ctypes.c_int,
          I64,I64,I64,I64,I64,ctypes.c_int,ctypes.c_int,I64]
    def array(q):
        L=[F(x)*FINE for x in q]
        assert all(x.denominator==1 for x in L)
        return np.array(list(map(int,L)),dtype=np.int64)
    root=array(tree.root["hi"])
    target=F(spec["target"])*2*Q*Q
    assert target.denominator>1
    def run(case,box,support):
        part,pos,d,n,old,new,depth,tried=case
        lo,hi=(array(z) for z in zip(*box,strict=True))
        out=np.zeros(29,dtype=np.int64)
        ret=guard.dynamic_guard_run(lo,hi,root,*T[d],ceil(target),depth+1,n,*support,
                                   len(support[1]),len(support[3])-1,out)
        assert ret==0 and out[:7].tolist()==[1,n,old,new,depth,0,tried],(part,pos,out[:7])
        assert n==2*(old+new)-1
        return case

    todo=[r for r in cases if not a.parts or r[0] in a.parts]
    assert todo
    result=[]
    with ThreadPoolExecutor(a.workers) as pool:
        queue=[]
        for q in todo:
            box=tree.leaf(*q[:2])
            queue.append(pool.submit(run,q,box,S.presolve(box)))
            if len(queue)>=a.workers*2:
                result.append(queue.pop(0).result())
        result.extend(f.result() for f in queue)
    assert result==todo
    print(json.dumps({"scope": "all 501 original retained leaves" if len(todo)==501 else "selected parts",
        "threshold":spec["target"],"input_sha256":digest(source),"localized_retained_sha256":spec["localization_retained_sha256"],
        "old_table_arrays_sha256":{n:r["canonical_arrays_sha256"] for n,r in spec["tables"].items()},
        "covered":len(result),"S_used_leaves":sum(q[5]>0 for q in result),"calls":sum(q[3] for q in result),
        "old_terminals":sum(q[4] for q in result),"S_terminals":sum(q[5] for q in result),"deepest":max(q[6] for q in result),
        "case_rows_sha256":digest(json.dumps(result,separators=(",",":")).encode()),
        "compiler":clang,"python":platform.python_version(),"python_flint":flint.__version__},indent=2))


if __name__=="__main__":
    main()
