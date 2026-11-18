"""Verify all B7 placements over a symmetric five-dimensional base box."""
from argparse import ArgumentParser
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from math import prod
from pathlib import Path
import json
import platform

import flint
import numpy as np
import certify_append_b7_at_three_minimum as frozen


def initialize(source,target,samples,depth,width):
    global limit,halfwidth
    frozen.initialize(source,target,samples,depth)
    limit,halfwidth=depth,width


def certify(job):
    num,box=job
    checked=0
    used=0
    unresolved=[]
    ledger=sha256()
    def visit(l,r,depth,path):
        nonlocal checked,used
        checked+=1
        c=frozen.FROZEN+tuple((a+b)/2 for a,b in zip(l,r))
        h=(halfwidth,)*5+tuple((b-a)/2 for a,b in zip(l,r))
        assert all(F(float(x))==x for x in c+h)
        c,h=np.array(c,dtype=float),np.array(h,dtype=float)
        passes,detail=frozen.checker.check_bound_leaf(c,h)
        if passes:
            used+=1
            ledger.update(json.dumps([path,detail]).encode()+b'\n')
            return
        if depth==limit:
            unresolved.append([list(map(str,l)),list(map(str,r))]);return
        axis=max(range(3),key=lambda j:(r[j]-l[j])*(F(1,2) if j==0 else 1))
        mid=(l[axis]+r[axis])/2
        visit(l,r[:axis]+(mid,)+r[axis+1:],depth+1,path+'0')
        visit(l[:axis]+(mid,)+l[axis+1:],r,depth+1,path+'1')
    visit(*box,0,'')
    assert checked==2*(used+len(unresolved))-1
    return {'index':num,'checked':checked,'taylor_checked':used,'residual':unresolved,'terminal_sha256':ledger.hexdigest()}


def main():
    parser=ArgumentParser();parser.add_argument('source',type=Path);parser.add_argument('output',type=Path)
    parser.add_argument('--width',type=F,default=F(1,512));parser.add_argument('--target',default='0.8363')
    parser.add_argument('--workers',type=int,default=2);parser.add_argument('--samples',type=int,default=512)
    parser.add_argument('--depth',type=int,default=16);parser.add_argument('--resume',type=Path);args=parser.parse_args()
    assert 0<args.width<F(1,64) and F(float(args.width))==args.width
    jobs=[]
    if args.resume:
        roots=json.loads(args.resume.read_bytes());assert roots['width']==str(args.width) and roots['target']==args.target
        jobs=[(j,tuple(tuple(map(F,r)) for r in box)) for j,box in enumerate(roots['residual'])]
    else:
        grid=(8,4,4);steps=[(b-a)/k for a,b,k in zip(frozen.ROOT_LO,frozen.ROOT_HI,grid)]
        for j,pos in enumerate(product(*(range(k) for k in grid))):
            l=tuple(a+k*h for a,k,h in zip(frozen.ROOT_LO,pos,steps));r=tuple(a+h for a,h in zip(l,steps))
            jobs.append((j,(l,r)))
        assert sum(prod(b-a for a,b in zip(*box)) for _,box in jobs)==prod(b-a for a,b in zip(frozen.ROOT_LO,frozen.ROOT_HI))
    rows=[];counts=Counter();remaining=[]
    with ProcessPoolExecutor(max_workers=args.workers,initializer=initialize,initargs=(args.source,args.target,args.samples,args.depth,args.width)) as pool:
        for j,row in enumerate(pool.map(certify,jobs,chunksize=1)):
            rows.append(row);counts['checked']+=row['checked'];counts['taylor']+=row['taylor_checked'];remaining.extend(row['residual'])
            if (j+1)%16==0:print(j+1,len(jobs),counts['checked'],len(remaining),flush=True)
    assert counts['checked']==2*(counts['taylor']+len(remaining))-len(jobs)
    report={'python':platform.python_version(),'flint':flint.__version__,'numpy':np.__version__,
        'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),'checker_sha256':frozen.CHECKER_SHA,
        'frozen_center':list(map(str,frozen.FROZEN)), 'width':str(args.width),'target':args.target,
        'samples':args.samples,'max_depth':args.depth,'jobs':len(jobs),'counts':dict(counts),
        'root':[list(map(str,frozen.ROOT_LO)),list(map(str,frozen.ROOT_HI))],
        'resume_sha256':sha256(args.resume.read_bytes()).hexdigest() if args.resume else None,
        'residual':remaining,'records':rows}
    args.output.write_text(json.dumps(report,separators=(',',':'))+'\n')
    print(args.target,str(args.width),dict(counts),len(remaining),flush=True)


if __name__=='__main__':
    main()
