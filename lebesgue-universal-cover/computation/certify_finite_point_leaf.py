"""Certify one complete Tomar five-body leaf after adjoining a fixed rational point set."""
import argparse
import ctypes
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
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
from audit_five_harmonic_tree import read_codes, split
from verify_tomar_fan_table import PINNED_MANIFEST, verified_table

FINE,Q,U=2**44,2**28,2**22
TABLE_SHA="fe5b816060d713642dc44e3f8981f4ae5d1b2b1c1e271b7dff75f031dd75d87c"


def ar(v):
    f=F(v)
    return arb(f.numerator)/f.denominator


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def original_leaf(tree,part,pos):
    manifest=tree/"SHA256SUMS"
    assert digest(manifest)==PINNED_MANIFEST
    pinned={n.removeprefix("./"):v for v,n in (r.split() for r in manifest.read_text().splitlines())}
    data=tree/f"part{part:04}.npy"
    hp=tree/f"part{part:04}.json"
    assert digest(data)==pinned[data.name] and digest(hp)==pinned[hp.name]
    header=json.loads(hp.read_bytes())
    assert header["part"]==part and header["split0"]==[2]*11 and header["ns"]==[3,5,7,9]
    box=list(zip(map(F,header["lo"]),map(F,header["hi"])))
    for k in range(11):
        box=split(box,k)[(part>>(10-k))&1]
    codes,prefix=read_codes(data.read_bytes())
    assert 0 <= pos < len(codes) and len(codes)==header["boxes"]
    i=0
    while i!=pos:
        assert codes[i]<11
        right=prefix.index(prefix[i],i+2)
        side=int(pos>=right)
        box=split(box,codes[i])[side]
        i=right if side else i+1
    assert codes[i]>=253
    return header,box,digest(data)


def upper_q(v):
    up=round(float((v*Q).mid()))+2
    assert ar(F(up,Q))>v
    return up


def old_tables(nums,box,pi):
    vectors=[(ar(F(x,U)),ar(F(y,U)),ar(F(x*x+y*y,U*U)).sqrt()) for x,y in nums]
    disk=[upper_q(length/2) for _,_,length in vectors]
    ends=[]
    for n,axis in [(3,None),(5,2),(7,5),(9,8)]:
        l,h=(F(0),F(0)) if axis is None else box[axis]
        halfcos=(pi/(2*n)).cos()
        radius=1/(2*halfcos)
        assert radius<ar(F(3,5))
        secant=Q*F(3,5)*(h-l)**2/8*F(U+2,U)
        error=-(-secant.numerator//secant.denominator)
        for angle in (l,h):
            vertices=[(ar(angle)+2*pi*k/n).sin_cos() for k in range(n)]
            for x,y,length in vectors:
                options=[]
                for s,c in vertices:
                    projection=x*c+y*s
                    options.append(radius*projection)
                    if not projection > -length*halfcos:
                        options.append(length+radius*projection)
                ends.append(max(upper_q(z) for z in options)+error)
    return disk,ends


def point_tables(data):
    points=data["points"]
    den=data["denominator"]
    assert type(den) is int and den>0 and all(len(p)==2 and all(type(x) is int for x in p) for p in points)
    cloud=set(map(tuple,points))
    assert len(cloud)==len(points) and {(x,-y) for x,y in cloud}==cloud
    dia2=max((a-x)**2+(b-y)**2 for (a,b),(x,y) in combinations(points,2))
    assert dia2<=den*den
    nums=data["directions"]
    assert data["direction_scale"]==U
    assert all(len(v)==2 and all(type(x) is int for x in v) and
               (U-2)**2<v[0]**2+v[1]**2<(U+2)**2 for v in nums)
    cells=data["cells"]
    assert len(cells)==384
    width=arb.pi()/len(cells)
    sinw,cosw=width.sin_cos()
    dual,pointers=[],[0]
    for j,triples in enumerate(cells):
        s,c=((2*j+1)*width).sin_cos()
        assert triples
        for row in triples:
            assert len(row)==6 and all(type(x) is int for x in row)
            ids,material=row[:3],row[3:]
            assert len(set(ids))==3 and all(0<=i<len(nums) for i in ids)
            assert all(0<=i<len(points) for i in material)
            v=[nums[i] for i in ids]
            w=[v[(k+1)%3][0]*v[(k+2)%3][1]-v[(k+1)%3][1]*v[(k+2)%3][0] for k in range(3)]
            if max(w)<0:
                w=[-t for t in w]
            assert min(w)>0 and all(sum(wt*z[d] for wt,z in zip(w,v))==0 for d in (0,1))
            total=sum(w)
            A,B=arb(0),arb(0)
            for wt,(a,b),mid in zip(w,v,material):
                x,y=points[mid]
                nx,ny=a*c+b*s,-a*s+b*c
                A+=wt*(nx*x+ny*y)
                B+=wt*(ny*x-nx*y)
            assert A>0
            bound=(A*cosw-abs(B)*sinw)/(total*den*U)
            low=int((bound*Q).lower().floor().unique_fmpz())
            assert ar(F(low,Q))<=bound
            dual.extend([*ids,*w,total,low])
        pointers.append(len(dual)//8)
    return dia2,pointers,dual


def finite_clouds():
    from certify_nine import arc_body, append_point, S
    assert S==10**12
    spec=[]
    den=S+10
    for n in (0,3,5,7,9):
        points=[]
        if n:
            arc_body(points,(648+n-1)//n,arb.pi(),n,arb(0),arb(0),arb(0))
        else:
            for i in range(1024):
                s,c=(2*arb.pi()*i/1024).sin_cos()
                append_point(points,c/2,s/2)
        pts=sorted(set(points))
        max_diameter=max((a-x)**2+(b-y)**2 for (a,b),(x,y) in combinations(pts,2))
        assert max_diameter<=den*den and {(x,-y) for x,y in pts}==set(pts)
        spec.append(dict(body=n,denominator=den,points=pts,diameter_squared_numerator=max_diameter))
    assert [len(row["points"]) for row in spec]==[1024,648,650,651,648]
    pi=F(22,7)
    assert arb.pi()<ar(pi)
    error=max(pi*pi/(8*648**2),pi*pi/(4*1024**2))+F(12,S)
    r=F(1,2)-error
    lower=(F(".83731")-pi*error*error)/(1+2*error/r)
    assert lower>F(".83730015182")
    return spec,error,lower


def run():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("cert","table","work","cache"):
        p.add_argument(k,type=Path)
    p.add_argument("--input",type=Path,default=HERE/"finite-point-leaf-input.json")
    p.add_argument("--finite-clouds",action="store_true",help="certify all six rational groups on the full leaf")
    a=p.parse_args()
    assert __debug__
    a.work.mkdir(exist_ok=False)
    a.cache.mkdir(parents=True,exist_ok=True)
    ctx.prec=160
    input_path=a.input
    data=json.loads(input_path.read_bytes())
    assert data["part"]==180 and data["pos"]==23955 and F(data["target"])==F(".83730")
    header,box,tree_sha=original_leaf(a.cert,data["part"],data["pos"])
    assert digest(a.table)==TABLE_SHA
    table,counts=verified_table(a.table,a.cert,header)
    # The certified five-body point is strictly inside this whole original leaf.
    from certify_five_body_contact_stationary import BOX, model, interval
    _,pose,_,_=model(*(interval(*z) for z in BOX))
    assert all(ar(l)<v and v<ar(h) for v,(l,h) in zip(pose,box,strict=True))
    ctx.prec=160
    diameter,pointers,dual=point_tables(data)
    disk,endpoints=old_tables(data["directions"],box,arb.pi())
    src=HERE/"finite_point_leaf_guard.cpp"
    binary=(a.work/"point-guard.dylib").resolve()
    compiler=Path("/usr/bin/clang++")
    flags=["-O3","-std=c++17","-shared","-fPIC","-fsanitize=signed-integer-overflow",
           "-fsanitize-undefined-trap-on-error","-fmodules-cache-path="+str(a.cache.resolve())]
    subprocess.run([compiler,*flags,src,"-o",binary],check=True)
    fun=ctypes.CDLL(str(binary)).dynamic_guard_run
    vec=lambda typ:np.ctypeslib.ndpointer(dtype=typ,ndim=1,flags="C_CONTIGUOUS")
    i64,i32=vec(np.int64),vec(np.int32)
    fun.argtypes=[i64,i64,i64,i32,i32,ctypes.c_int64,ctypes.c_int,ctypes.c_int,
                  i64,i64,i64,i64,i64,ctypes.c_int,ctypes.c_int,i64]
    fun.restype=ctypes.c_int
    ints=lambda seq:np.ascontiguousarray(seq,dtype=np.int64).ravel()
    params=[[FINE*row[j] for row in box] for j in (0,1)]+[[FINE*F(x) for x in header["hi"]]]
    assert all(F(x).denominator==1 for r in params for x in r)
    target=data["target"]
    expected=data["expected"]
    if a.finite_clouds:
        assert digest(a.input)==digest(HERE/"finite-point-leaf-input.json")
        target=".83731"
        expected=[1,3182625,47209,1544104,83,0,3125697]
    threshold=F(target)*2*Q*Q
    goal=threshold.numerator//threshold.denominator+1
    stats=np.zeros(29,np.int64)
    assert fun(*map(ints,params),table.ravel(),counts,goal,data["depth_limit"],data["max_calls"],
               ints(data["directions"]),ints(disk),ints(endpoints),ints(pointers),ints(dual),
               len(disk),len(pointers)-1,stats)==0
    ok,calls,old,new,depth,reason,attempts=map(int,stats[:7])
    assert [ok,calls,old,new,depth,reason,attempts]==expected
    assert ok==1 and calls==2*(old+new)-1 and reason==0 and new>0
    out=dict(part=data["part"],pos=data["pos"],target=target,goal_integer=goal,
             points=len(data["points"]),denominator=data["denominator"],diameter_squared_numerator=diameter,
             angular_cells=len(pointers)-1,balanced_triples=len(dual)//8,verified_reflection=True,
             calls=calls,old_terminals=old,point_terminals=new,depth=depth,certified_root_contains_contact=True,
             pinned_source_manifest=PINNED_MANIFEST,original_tree_sha256=tree_sha,
             table_sha256=TABLE_SHA,driver_sha256=digest(Path(__file__)),cpp_sha256=digest(src),
             input_sha256=digest(input_path),binary_sha256=digest(binary),compiler_flags=flags,
             python=platform.python_version(),flint=flint.__version__)
    if a.finite_clouds:
        groups,error,lower=finite_clouds()
        groups.append(dict(body="P",denominator=data["denominator"],points=data["points"],
                           diameter_squared_numerator=diameter))
        cloud_path=a.work/"six-rational-groups.json"
        cloud_path.write_text(json.dumps(groups,separators=(",",":"))+"\n")
        out.update(finite_group_sizes=[len(g["points"]) for g in groups],cloud_sha256=digest(cloud_path),
                   hausdorff_error=str(error),all_finite_area_lower=str(lower),all_finite_decimal_floor=".83730015182")
    (a.work/"replay.json").write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps(out,indent=2))


if __name__=="__main__":
    run()
