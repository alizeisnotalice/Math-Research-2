"""Compare P30, the original Q11 guard, and certified mixed-area center guards."""
from argparse import ArgumentParser
import ctypes
from fractions import Fraction as F
import json
from pathlib import Path
import subprocess
import sys
from time import process_time

import numpy as np
from flint import arb, ctx


def intersections(guards,tree,index):
    levels=read(index)
    hits=[]
    for part,positions in levels.items():
        for pos in positions:
            ends=[(int(a*2**45),int(b*2**45)) for a,b in tree.leaf(int(part),pos)]
            for row in guards:
                d=[max(a-z,z-b,0) for (a,b),z in zip(ends,row[:-1],strict=True)]
                fits=True
                for j,l in ((0,0),(3,163),(6,115),(9,89)):
                    room=1000*row[-1]-(l*d[j-1] if j else 0)
                    if room<0 or 10**6*(d[j]**2+d[j+1]**2)>room*room:
                        fits=False
                        break
                if fits:
                    hits.append([int(part),pos])
                    break
    return dict(overlap=len(hits),all_leaves=sum(map(len,levels.values())),leaves=hits)


def main(a):
    ctx.prec=160
    sys.path.insert(0,str(a.comp))
    from build_q11_leaf_bank import reflected_average, EPS
    from replay_q11_with_p30 import PointProof
    from certify_tomar_six_high import FINE, Q, U, SourceTree, build_kernel, ceil, digest, load_table
    from certify_nine import ar

    source=a.comp/'tomar_six_083730_guard.cpp'
    assert digest(source.read_bytes())=="42a4ee583d22940b3cd0b05f04285ba7e3d351745de707476211dce4efba5f81"
    pi=arb.pi()
    r=ar(F(".707"))
    assert (pi-(1/(2*r)).acos())/4+(r*r-ar(F(1,4))).sqrt()/2 > ar(a.target)
    assert r*r < ar(F(1,2))
    data=read(a.comp/'tomar-six-high-input.json')
    tree=SourceTree(a.tree)
    table=load_table(a.table,tree,data['tables']['64'])
    pcloud=read(a.comp/'finite-point-leaf-input.json')
    pcloud['cells']=read(a.comp/'finite-point-global-input.json')['cells']
    pp=PointProof(pcloud)
    qcloud,qden=reflected_average(read(a.comp/'multi-host-points.json'))
    bank=read(a.comp/'q11-leaf-bank.json')
    qq=PointProof({**bank,'points':qcloud,'denominator':qden,'direction_scale':U},EPS)
    for n,c in ((5,163),(7,115),(9,89)):
        assert (pi/(2*n)).tan()/2 < ar(F(c,1000))
    work=a.work.resolve()
    build_kernel(work,data,a.cache)
    lib=work/'guard.so'
    extra=Path(__file__).with_name('mixed_point_leaf_guard.cpp')
    subprocess.run(['clang++','-O3','-std=c++17','-ftrapv','-fPIC','-shared',
                    f'-fmodules-cache-path={a.cache}','-I',str(work),'-I',str(a.comp),
                    str(extra),'-o',str(lib)],check=True)
    call=ctypes.CDLL(str(lib)).mixed_point_leaf_run
    vec=lambda d:np.ctypeslib.ndpointer(dtype=d,ndim=1,flags='C_CONTIGUOUS')
    i64,i32,M=vec(np.int64),vec(np.int32),ctypes.c_int
    call.argtypes=[i64,i64,i64,i32,i32,ctypes.c_int64,M,M,M,*[i64]*5,M,*[i64]*5,M,i64,M,i64,M,M,M,i64]
    call.restype=M

    target=ceil(a.target*2*Q*Q)
    assert F(target,2*Q*Q)>a.target
    def scaled(vals):
        rationals=[F(z)*FINE for z in vals]
        assert all(z.denominator==1 for z in rationals)
        return np.array(list(map(int,rationals)),dtype=np.int64)
    original=scaled(tree.root['hi'])
    before=process_time()
    if a.cert is not None:
        proposed=read(a.cert)
        if proposed.get('angle_verifier')=='events':
            sys.path.insert(0,str(Path(__file__).resolve().parent))
            from compress_q11_event_banks import replay
            certs=replay(a.comp,proposed)
        else:
            from certify_mixed_leaf_centers import replay_bank
            certs=replay_bank(a.comp,proposed)
    else:
        certs=None
    bank_cpu=process_time()-before
    if certs:
        assert F(certs['target'])==a.target and certs['scale']==2*FINE
        assert certs['Q11_sha256']==digest((a.comp/'multi-host-points.json').read_bytes())
    if a.overlap_index is not None:
        assert certs and a.reuse
        print(json.dumps(intersections(certs['leaves'][a.reuse]['guards'],tree,a.overlap_index)),flush=True)
    for key in a.cases:
        box=tree.leaf(*map(int,key.split('/')))
        lo,hi=(scaled(t) for t in zip(*box,strict=True))
        before=process_time()
        d1,d2=pp.presolve(box),qq.presolve(box)
        prep_cpu=process_time()-before
        records=certs['leaves'].get(a.reuse or key,{}) if certs else {}
        material=[r+[int(c.get('strategy','mixed')=='cap')]
                  for r,c in zip(records.get('guards',[]),records.get('certs',[]),strict=True)]
        coords=np.array(material,np.int64).ravel()
        assert len(coords)%13==0
        for mode in a.modes:
            if mode&2:
                assert len(coords)>0
            out=np.zeros(35+len(material),np.int64)
            trace=np.zeros(22*a.trace,np.int64)
            before=process_time()
            assert call(lo,hi,original,*table,target,a.depth,a.budget,mode,
                        *d1,len(d1[1]),*d2,len(d2[1]),coords,len(coords)//13,
                        trace,a.trace,a.stride,a.trace_kind,out)==0
            used=process_time()-before
            stats=out[:35].tolist()
            uses=out[35:].tolist()
            assert stats[5]==sum(stats[33:35])==sum(uses)
            if stats[0]:
                assert stats[1]==2*sum(stats[2:6])-1 and stats[7]==0
            if not mode&1:
                assert stats[4]==0
            if not mode&2:
                assert stats[5]==0
            row=dict(leaf=key,goal=str(a.target),mode=mode,budget=a.budget,stats=stats,
                     prep_cpu=prep_cpu,search_cpu=used,bank_cpu=bank_cpu,guards=len(coords)//13,guard_uses=uses)
            print(json.dumps(row),flush=True)
            if a.trace and (mode&1 or a.trace_kind==1):
                a.work.joinpath(key.replace('/','-')+f'-{mode}-trace.json').write_text(
                    json.dumps(trace[:22*stats[10]].reshape(-1,22).tolist(),separators=(',',':'))+'\n')


def read(path):
    return json.loads(path.read_bytes())


if __name__=='__main__':
    p=ArgumentParser(description=__doc__)
    for key in ('comp','tree','table','work','cache'):
        p.add_argument(key,type=Path)
    p.add_argument('--cert',type=Path)
    p.add_argument('--overlap-index',type=Path)
    p.add_argument('--reuse',help='reuse another leaf\'s certified geometric guards')
    p.add_argument('--target',type=F,required=True)
    p.add_argument('--budget',type=int,required=True)
    p.add_argument('--depth',type=int,default=116)
    p.add_argument('--modes',type=int,choices=(0,1,2,3),nargs='+',required=True)
    p.add_argument('--trace',type=int,default=0)
    p.add_argument('--trace-kind',type=int,choices=(0,1),default=0)
    p.add_argument('--stride',type=int,default=1)
    p.add_argument('--cases',nargs='+',required=True)
    args=p.parse_args()
    assert args.stride>0 and args.trace>=0
    main(args)
