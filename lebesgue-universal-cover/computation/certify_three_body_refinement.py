"""Arb replay of a tighter three-body threshold on exact clipped dyadic seed cells."""
from argparse import ArgumentParser
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from io import BytesIO
from itertools import product
import json
from pathlib import Path
import platform
import tarfile

import flint
from flint import arb
import numpy as np

SOURCE_COMMIT='10af98606c98ca9ffd0d9f402cae1e2efe74b792'
CHECKER_SHA='a53e3b93c87dfb6bec28b4476989e008c59f4d28ec0bdbe3dc284f9beeb39051'
SEED_SHA='2026ebe7e5339b89542e0429e5ad704b1a956278e3ae455b6d667b4c11892846'
ROOT_LO=tuple(map(Q,[-9/128,-9/128,0,-3/32,-3/32]))
ROOT_HI=tuple(map(Q,[9/128,9/128,27/128,3/32,3/32]))


def split(lo,hi,axis):
    m=(lo[axis]+hi[axis])/2
    return (lo,hi[:axis]+(m,)+hi[axis+1:]),(lo[:axis]+(m,)+lo[axis+1:],hi)


def determinant(a,b):
    return a[0]*b[1]-a[1]*b[0]


def corner_bounds(runs,pa,h):
    """Bound one verified fan on the full four-translation/one-rotation cell."""
    def twice_area(cs,sn,translations):
        def rot(body,p):
            if body!=2:
                return p
            return cs*p[0]-sn*p[1],sn*p[0]+cs*p[1]
        z=arb(0)
        for i,(body,start,end,intrinsic) in enumerate(runs):
            start,end=rot(body,start),rot(body,end)
            t=translations[body]; b2,next_start=runs[(i+1)%len(runs)][:2]
            u=translations[b2]; v=rot(b2,next_start)
            z+=intrinsic+determinant(t,(end[0]-start[0],end[1]-start[1]))
            z+=determinant((t[0]+end[0],t[1]+end[1]),(u[0]+v[0],u[1]+v[1]))
        return z
    mid,rad=pa[2][:2]; cm,sm=mid.cos(),mid.sin()
    for signs in product([-1,1],repeat=4):
        t=[(arb(0),arb(0))]
        for body,offset in ((1,0),(2,3)):
            center=pa[body][2]; first=0 if body==1 else 2
            t.append(tuple(center[k]+arb(float(h[offset+k]))*signs[first+k] for k in range(2)))
        f0=twice_area(1,0,t); fpi=twice_area(-1,0,t)
        a=(f0+fpi)/2; b=(f0-fpi)/2; c=twice_area(0,1,t)-a
        lower=a+b*cm+c*sm-abs(-b*sm+c*cm)*rad-(abs(b)+abs(c))*rad**2/2
        yield lower/2


def corner_ok(runs,pa,h,tau):
    return all(lower>tau for lower in corner_bounds(runs,pa,h))


def part_run(job):
    source,part,target,max_extra=job; source=Path(source)
    verifier=source/'checker/lbcheck.py'
    assert sha256(verifier.read_bytes()).hexdigest()==CHECKER_SHA
    spec=spec_from_file_location('tomar_lbcheck',verifier); lib=module_from_spec(spec); spec.loader.exec_module(lib)
    seed=source/'certificates/cert5d_08347.tar.xz'; assert sha256(seed.read_bytes()).hexdigest()==SEED_SHA
    with tarfile.open(seed,'r:xz') as t:
        def read(ext): return t.extractfile(f'cert5d_08347/part{part:04d}.{ext}').read()
        meta=json.loads(read('json')); codes=np.load(BytesIO(read('npy')), allow_pickle=False)
    assert meta['ns']==[3,5] and meta['split0']==[4]*5 and meta['part']==part and meta['nparts']==16
    lo,hi=tuple(map(Q,meta['lo'])),tuple(map(Q,meta['hi']))
    assert all(a<=b<c<=d for a,b,c,d in zip(lo,ROOT_LO,ROOT_HI,hi))
    tau=arb(target); counts=Counter(); covered=Q(0); unresolved=[]; deepest=0
    checker=lib.Checker([3,5],float(target),meta['N'],meta['eps']); checker.tau=tau
    original=checker.taylor_ok
    def test(runs,pa,h):
        valid,msg=original(runs,pa,h)
        if valid:
            checker.method='taylor'; return True,msg
        if corner_ok(runs,pa,h,tau):
            checker.method='corners'; return True,'corners'
        return False,'both lower bounds inconclusive'
    checker.taylor_ok=test
    def examine(l,r,depth):
        nonlocal deepest
        deepest=max(deepest,depth); counts['checked']+=1
        c=[(a+b)/2 for a,b in zip(l,r)]; w=[(b-a)/2 for a,b in zip(l,r)]
        assert all(Q(float(x))==x for x in c+w)
        valid,_=checker.check_bound_leaf(np.array(c,dtype=float),np.array(w,dtype=float))
        if valid:
            counts[checker.method]+=1; return
        if depth==max_extra:
            unresolved.append([[str(a) for a in l],[str(b) for b in r]]);return
        # Only the choice of splitting coordinate is heuristic.
        axis=max(range(5),key=lambda i:float(w[i])*(0.22 if i==2 else 1))
        for ll,rr in split(l,r,axis): examine(ll,rr,depth+1)
    pos=0
    def visit(l,r):
        nonlocal pos,covered
        code=int(codes[pos]);pos+=1
        if 0<=code<5:
            for ll,rr in split(l,r,code): visit(ll,rr)
            return
        assert code in (254,255)
        l=tuple(map(max,l,ROOT_LO)); r=tuple(map(min,r,ROOT_HI))
        if any(a>=b for a,b in zip(l,r)): return
        # The closures of positive-volume clipped tiles cover the closed root.
        v=Q(1)
        for a,b in zip(l,r):v*=b-a
        covered+=v; counts['seed_tiles']+=1;examine(l,r,0)
    for grid in list(product(range(4),repeat=5))[part::16]:
        sz=[(b-a)/4 for a,b in zip(lo,hi)]
        visit(tuple(a+j*z for a,j,z in zip(lo,grid,sz)),tuple(a+(j+1)*z for a,j,z in zip(lo,grid,sz)))
    assert pos==len(codes)
    return {'part':part,'seed_root':[[str(x) for x in lo],[str(x) for x in hi]],'volume':str(covered),
            'counts':dict(counts),'extra_depth':deepest,'unresolved':unresolved}


def main():
    parser=ArgumentParser();parser.add_argument('source');parser.add_argument('--workers',type=int,default=3)
    parser.add_argument('--threshold',default='0.83475');parser.add_argument('--max-extra',type=int,default=15)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    assert __debug__ and 0<Q(args.threshold)<=Q(419,500)
    assert Q(7,100)<ROOT_HI[0] and Q(227,2500)<ROOT_HI[3] and Q(22,7)/15<ROOT_HI[2]
    jobs=[(args.source,i,args.threshold,args.max_extra) for i in range(16)]
    with ProcessPoolExecutor(max_workers=args.workers) as ex: parts=list(ex.map(part_run,jobs))
    assert len({json.dumps(p['seed_root']) for p in parts})==1
    volume=Q(1)
    for a,b in zip(ROOT_LO,ROOT_HI): volume*=b-a
    assert sum(Q(p['volume']) for p in parts)==volume
    counts=sum((Counter(p['counts']) for p in parts),Counter());unresolved=sum(len(p['unresolved']) for p in parts)
    status={'target':args.threshold,'completed':unresolved==0, 'unresolved':unresolved,
            'counts':dict(counts),'root':[[str(x) for x in r] for r in [ROOT_LO,ROOT_HI]],'volume':str(volume),
            'python':platform.python_version(),'numpy':np.__version__,'flint':flint.__version__,
            'source_commit':SOURCE_COMMIT,'checker_sha256':CHECKER_SHA,'seed_sha256':SEED_SHA,
            'max_extra':args.max_extra,'parts':parts}
    args.output.write_text(json.dumps(status,indent=2)+'\n'); print(json.dumps({k:status[k] for k in ['completed','target','counts','unresolved']}),flush=True)
    # The second-stage checker handles every explicitly archived residual cell.

if __name__=='__main__': main()
