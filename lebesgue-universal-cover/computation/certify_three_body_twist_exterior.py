from argparse import ArgumentParser
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as Q
from hashlib import sha256
from io import BytesIO
from itertools import product
from math import prod
import json
from pathlib import Path
import sys
import tarfile
import numpy as np

parser=ArgumentParser()
parser.add_argument('source',type=Path)
parser.add_argument('band',choices=('wide','near'))
parser.add_argument('output',type=Path)
parser.add_argument('--computation',type=Path,default=Path(__file__).parent)
parser.add_argument('--workers',type=int,default=2)
args=parser.parse_args()
sys.path.insert(0,str(args.computation))
import certify_three_body_weighted as weighted
from certify_three_body_refinement import CHECKER_SHA,ROOT_LO,ROOT_HI,SEED_SHA,split

source=args.source
target='0.834780945776'
D=65536
if args.band=='wide':
    angular_lo,angular_hi=Q(0),Q(163,D)
    inner_lo=(-1114,-294,0,1213,-294)
    inner_hi=(-499,294,163,1795,294)
else:
    angular_lo,angular_hi=Q(163,D),Q(294,D)
    inner_lo=(-950,-49,163,1436,-100)
    inner_hi=(-617,49,294,1537,25)
front_lo=(*ROOT_LO[:2],angular_lo,*ROOT_LO[3:])
front_hi=(*ROOT_HI[:2],angular_hi,*ROOT_HI[3:])
remove_lo=tuple(Q(v,D) for v in inner_lo)
remove_hi=tuple(Q(v,D) for v in inner_hi)


def intersect(l,r,a,b):
    return tuple(map(max,l,a)),tuple(map(min,r,b))


def volume(l,r):
    values=[b-a for a,b in zip(l,r)]
    return prod(values) if min(values)>0 else Q(0)


def subtract(l,r,a,b):
    if volume(*intersect(l,r,a,b))==0:
        yield (l,r)
        return
    for k in (0,1,3,4):
        if l[k]<a[k]:
            yield l,r[:k]+(a[k],)+r[k+1:]
            l=l[:k]+(a[k],)+l[k+1:]
        if b[k]<r[k]:
            yield l[:k]+(b[k],)+l[k+1:],r
            r=r[:k]+(b[k],)+r[k+1:]


def seed_jobs():
    seed=source/'certificates/cert5d_08347.tar.xz'
    assert sha256(seed.read_bytes()).hexdigest()==SEED_SHA
    jobs=[]
    with tarfile.open(seed,'r:xz') as t:
        for part in range(16):
            def read(ext):return t.extractfile(f'cert5d_08347/part{part:04d}.{ext}').read()
            meta=json.loads(read('json'));a,b=tuple(map(Q,meta['lo'])),tuple(map(Q,meta['hi']));sz=[(y-x)/4 for x,y in zip(a,b)]
            codes=np.load(BytesIO(read('npy')),allow_pickle=False)
            assert meta['ns']==[3,5] and meta['split0']==[4]*5 and meta['part']==part and meta['nparts']==16
            assert all(v<=x<y<=w for v,x,y,w in zip(a,front_lo,front_hi,b))
            pos=0
            def visit(l,r):
                nonlocal pos
                code=int(codes[pos]);pos+=1
                if 0<=code<5:
                    for ll,rr in split(l,r,code):visit(ll,rr)
                    return
                assert code in (254,255)
                l,r=intersect(l,r,front_lo,front_hi)
                if volume(l,r)==0:return
                for ll,rr in subtract(l,r,remove_lo,remove_hi):
                    assert volume(ll,rr)>0 and volume(*intersect(ll,rr,remove_lo,remove_hi))==0
                    jobs.append((part,len(jobs),(ll,rr)))
            for grid in list(product(range(4),repeat=5))[part::16]:
                visit(tuple(x+j*d for x,j,d in zip(a,grid,sz)),tuple(x+(j+1)*d for x,j,d in zip(a,grid,sz)))
            assert pos==len(codes)
    assert sum(volume(*box) for _,_,box in jobs)+volume(remove_lo,remove_hi)==volume(front_lo,front_hi)
    return jobs


def work(job):
    _,_,(l,r)=job
    cf,hf=[(a+b)/2 for a,b in zip(l,r)],[(b-a)/2 for a,b in zip(l,r)]
    assert all(Q(float(a))==a for a in cf+hf)
    c,h=np.array(cf,dtype=float),np.array(hf,dtype=float)
    if weighted.checker.check_bound_leaf(c,h)[0]:
        return {'counts':{'checked':1,'taylor':1},'unresolved':[]}
    return weighted.solve(job)


def main():
    jobs=seed_jobs(); print('jobs',len(jobs),flush=True)
    counts=Counter();bad=[];mx=0;ld=None
    ledger=sha256()
    with ProcessPoolExecutor(max_workers=args.workers,initializer=weighted.init_worker,initargs=(source,target,10)) as ex:
        for i,row in enumerate(ex.map(work,jobs,chunksize=1)):
            counts.update(row['counts']);mx=max(mx,row.get('depth',0))
            if row.get('min_bound') is not None:
                v=Q(row['min_bound']);ld=v if ld is None else min(ld,v)
            bad.extend(row['unresolved']);ledger.update(json.dumps(row,sort_keys=True).encode()+b'\n')
            if (i+1)%50==0:print(i+1,counts['checked'],len(bad),flush=True)
    tally=sum(v for k,v in counts.items() if k=='taylor' or k.startswith('weights_'))+len(bad)
    assert counts['checked']==2*tally-len(jobs)
    data={'target':target,'band':args.band,'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
          'seed_sha256':SEED_SHA,'checker_sha256':CHECKER_SHA,'ledger_sha256':ledger.hexdigest(),'frontier':[list(map(str,front_lo)),list(map(str,front_hi))],'excluded':[list(map(str,remove_lo)),list(map(str,remove_hi))],
          'jobs':len(jobs),'counts':dict(counts),'deepest':mx,'weighted_min':str(ld),'residual':bad}
    args.output.write_text(json.dumps(data,indent=2)+'\n');print({k:v for k,v in data.items() if k!='residual'},flush=True)

if __name__=='__main__':main()
