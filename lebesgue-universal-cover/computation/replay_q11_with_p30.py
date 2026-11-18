"""Compare P30 with sequential independent P30/Q11 exclusions on exact original boxes."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import ctypes
from fractions import Fraction as F
import json
from pathlib import Path
import subprocess

import numpy as np
from flint import arb, ctx

from build_q11_leaf_bank import reflected_average, EPS
from certify_finite_point_leaf import point_tables, upper_q
from certify_nine import ar
from certify_tomar_six_high import FINE, Q, U, SourceTree, DirectionProof, build_kernel, ceil, digest, load_table

HERE=Path(__file__).resolve().parent

class PointProof(DirectionProof):
    def __init__(self,d,error=F(0)):
        self.pi=arb.pi()
        diameter,pointers,rows=point_tables(d)
        self.v=[(ar(F(x,U)),ar(F(y,U)),ar(F(x*x+y*y,U*U)).sqrt()) for x,y in d["directions"]]
        self.directions=np.asarray(d["directions"],dtype=np.int64).ravel()
        self.disk=np.asarray([upper_q(v[2]/2) for v in self.v],dtype=np.int64)
        self.pointers=np.asarray(pointers,dtype=np.int64)
        self.triples=np.asarray(rows,dtype=np.int64)
        self.triples[7::8]-=ceil(Q*F(U+2,U)*error)

def read(path):
    return json.loads(path.read_bytes())

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('comp','tree','table','bank','work','cache'):
        p.add_argument(name,type=Path)
    p.add_argument('--cases',nargs='+')
    p.add_argument('--bits',type=int,default=64)
    p.add_argument('--budget',type=int,default=1048575)
    p.add_argument('--depth',type=int,default=116)
    p.add_argument('--workers',type=int,default=2)
    p.add_argument('--mode',choices=('both','P30','Q11'),default='both')
    a=p.parse_args()
    ctx.prec=160
    pi=arb.pi()
    r=ar(F(".707"))
    assert (pi-(1/(2*r)).acos())/4+(r*r-ar(F(1,4))).sqrt()/2 > ar(F(".838"))
    assert r*r < ar(F(1,2))
    data=read(a.comp/'tomar-six-high-input.json')
    tree=SourceTree(a.tree)
    table=load_table(a.table,tree,data['tables'][str(a.bits)])
    pcloud=read(a.comp/'finite-point-leaf-input.json')
    pcloud["cells"]=read(a.comp/'finite-point-global-input.json')['cells']
    pp=PointProof(pcloud)
    bank=read(a.bank)
    qcloud,qden=reflected_average(read(a.comp/'multi-host-points.json'))
    assert len(qcloud)==11
    qq=PointProof({**bank,"points":qcloud,"denominator":qden,"direction_scale":U},EPS)

    work,cache=a.work.resolve(),a.cache.resolve()
    build_kernel(work,data,cache)
    mode2=a.comp/'tomar_six_083730_guard.cpp'
    assert digest(mode2.read_bytes())=="42a4ee583d22940b3cd0b05f04285ba7e3d351745de707476211dce4efba5f81"
    lib=work/'p30-q11.so'
    subprocess.run(['clang++','-O3','-std=c++17','-ftrapv','-fPIC','-shared',
          f'-fmodules-cache-path={cache}','-I',str(work),'-I',str(a.comp),
          str(HERE/'dual_point_guard.cpp'),'-o',str(lib)],check=True)
    guard=ctypes.CDLL(str(lib)).dual_point_guard_run
    vect=lambda d:np.ctypeslib.ndpointer(dtype=d,ndim=1,flags='C_CONTIGUOUS')
    i64,i32,M=vect(np.int64),vect(np.int32),ctypes.c_int
    guard.argtypes=[i64,i64,i64,i32,i32,ctypes.c_int64,M,M,M,*[i64]*5,M,*[i64]*5,M,i64]
    guard.restype=M
    caps=ceil(F('.83740')*2*Q*Q)
    assert F(caps,2*Q*Q)>F('.83740')
    def scaled(vals):
        rationals=[F(z)*FINE for z in vals]
        assert all(z.denominator==1 for z in rationals)
        return np.array(list(map(int,rationals)),dtype=np.int64)

    source=scaled(tree.root['hi'])
    keys=[tuple(map(int,c.split('/'))) for c in a.cases] if a.cases else list(map(tuple,bank['cases']))
    modes={'both':[0,1],'P30':[0],'Q11':[1]}[a.mode]

    def run(key,mode,lo,hi,d1,d2):
        ans=np.zeros(33,dtype=np.int64)
        assert 0==guard(lo,hi,source,*table,caps,a.depth,a.budget,mode,
            *d1,len(d1[1]),*d2,len(d2[1]),ans)
        z=ans.tolist()
        if z[0]:
            assert z[1]==2*sum(z[2:5])-1 and not z[6]
        if not mode:
            assert z[4]==z[8]==z[10]==0
        row={'part':key[0],'pos':key[1],'table':a.bits,'budget':a.budget,'depth_limit':a.depth,
             'target':'.83740','Q11_enabled':bool(mode),'stats':z}
        print(json.dumps(row),flush=True)

    with ThreadPoolExecutor(a.workers) as pool:
        futures=[]
        for key in keys:
            box=tree.leaf(*key)
            lo,hi=(scaled(z) for z in zip(*box,strict=True))
            d1,d2=pp.presolve(box),qq.presolve(box)
            for mode in modes:
                futures.append(pool.submit(run,key,mode,lo,hi,d1,d2))
        for result in futures:
            result.result()

if __name__=='__main__':
    main()
