"""Three-dimensional inner-fan proof at the (rounded) exact three-body minimum."""
from argparse import ArgumentParser
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec,spec_from_file_location
from itertools import product
import json
from math import floor,prod
from pathlib import Path
import platform
import warnings

import flint
from flint import arb
import numpy as np
from scipy.optimize import OptimizeWarning,linprog

SCALE=10**12
CHECKER_SHA='a53e3b93c87dfb6bec28b4476989e008c59f4d28ec0bdbe3dc284f9beeb39051'
DEN=2**42
A=F(round(F('-.0121545587154054085')*DEN),DEN)
B=F(round(F('.02299799824774089')*DEN),DEN)
FROZEN=(A,F(0),F(0),B,F(0))
ROOT_LO=(F(0),F(-1,16),F(-1,16))
ROOT_HI=(F(29,64),F(1,16),F(1,16))
SIGNS=list(product((-1,1),repeat=2))


def determinant(a,b):
    return a[0]*b[1]-a[1]*b[0]


def row_lower(runs,pa,h):
    def twice_area(cs,sn,trans):
        def pt(owner,p):
            if owner==3:
                x,y=p;p=(cs*x-sn*y,sn*x+cs*y)
            t=trans[owner]
            return t[0]+p[0],t[1]+p[1]
        value=arb(0)
        for (i,x,y,z),(j,w,_,_) in zip(runs,runs[1:]+runs[:1]):
            # Each within-body contribution has a determinant invariant under rotation.
            x,y=pt(i,x),pt(i,y); t=trans[i]
            value+=z+determinant(t,(y[0]-x[0],y[1]-x[1]))+determinant(y,pt(j,w))
        return value
    theta,dtheta=pa[3][:2]
    cs,sn=theta.cos(),theta.sin()
    initial=[pa[i][2] for i in range(3)]
    ret=[]
    for sx,sy in SIGNS:
        trans=initial+[(pa[3][2][0]+arb(h[6])*sx,pa[3][2][1]+arb(h[7])*sy)]
        z0,zpi=twice_area(1,0,trans),twice_area(-1,0,trans)
        a,b=(z0+zpi)/2,(z0-zpi)/2
        c=twice_area(0,1,trans)-a
        low=(a+b*cs+c*sn-abs(-b*sn+c*cs)*dtheta-(abs(b)+abs(c))*dtheta**2/2)/2
        value=floor(float(low.lower())*SCALE)-1
        assert low>arb(value)/SCALE
        ret.append(value)
    return ret


def weight_ok(rows,weights):
    total=sum(weights)
    return total>0 and all(w>=0 for w in weights) and all(
        sum(w*row[j] for w,row in zip(weights,rows))>target_integer*total for j in range(4))


def dual(rows):
    for i,row in enumerate(rows):
        if min(row)>target_integer:
            return [1 if j==i else 0 for j in range(len(rows))]
    if not rows or any(max(row[v] for row in rows)<=target_integer for v in range(4)):
        return None
    last=rows[-1]
    for j,row in enumerate(rows[:-1]):
        lo,hi=F(0),F(1)
        for x,y in zip(row,last):
            if x>y:
                lo=max(lo,F(target_integer-y,x-y))
            elif x<y:
                hi=min(hi,F(target_integer-y,x-y))
            elif x<=target_integer:
                break
        else:
            if lo<hi:
                w=(lo+hi)/2
                weights=[0]*len(rows)
                weights[j]=w.numerator
                weights[-1]=w.denominator-w.numerator
                if weight_ok(rows,weights):
                    return weights
    if len(rows)<3:
        return None
    gap=(np.array(rows,dtype=np.int64)-target_integer)/1e6
    n=len(rows)
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore',category=OptimizeWarning,message='Unrecognized options detected.*')
        result=linprog(np.r_[np.zeros(n),-1.],A_ub=np.column_stack([-gap.T,np.ones(4)]),b_ub=np.zeros(4),
            A_eq=[np.r_[np.ones(n),0.]],b_eq=[1.],bounds=[(0,None)]*n+[(None,None)],method='highs',options={'threads':1})
    if not result.success:
        return None
    weights=[round(float(w)*10**8) for w in result.x[:-1]]
    return weights if weight_ok(rows,weights) else None


def initialize(source,target,samples,depth):
    global checker,target_integer,depth_limit
    file=source/'checker/lbcheck.py'
    assert sha256(file.read_bytes()).hexdigest()==CHECKER_SHA
    spec=spec_from_file_location('tomar_lbcheck',file);module=module_from_spec(spec);spec.loader.exec_module(module)
    checker=module.Checker([3,5,7],float(target),samples,3e-6)
    checker.tau=arb(target)
    ti=F(target)*SCALE
    assert ti.denominator==1
    target_integer=ti.numerator
    depth_limit=depth


def process(job):
    index,box=job
    stats=Counter(); residual=[]; ledger=sha256(); lower=None
    def visit(l,r,depth,trail):
        nonlocal lower
        stats['checked']+=1
        cq=FROZEN+tuple((a+b)/2 for a,b in zip(l,r));hq=(F(0),)*5+tuple((b-a)/2 for a,b in zip(l,r))
        assert all(F(float(x))==x for x in cq+hq)
        c,h=np.array(cq,dtype=float),np.array(hq,dtype=float)
        if checker.check_bound_leaf(c,h)[0]:
            stats['taylor']+=1;ledger.update((trail+'t\n').encode());return
        rows=[]
        attempts=[(0,0,0),(-1,0,0),(1,0,0),*((0,sx,sy) for sx,sy in SIGNS)]
        for shift in attempts:
            sel=np.array([0]*5+list(shift));runs,pa=checker.verified_runs(c,h,c+sel*h)
            if runs is None: continue
            rows.append(row_lower(runs,pa,h))
            weights=dual(rows)
            if weights is None: continue
            assert weight_ok(rows,weights)
            stats['weights_'+str(sum(w>0 for w in weights))]+=1
            val=min(F(sum(w*row[j] for w,row in zip(weights,rows)),sum(weights)*SCALE) for j in range(4))
            assert val>F(target_integer,SCALE)
            lower=val if lower is None else min(lower,val)
            ledger.update(json.dumps([trail,rows,weights,str(val)]).encode()+b'\n');return
        if depth>=depth_limit:
            residual.append([list(map(str,l)),list(map(str,r))]);return
        axis=max(range(3),key=lambda j:(r[j]-l[j])*(F(1,2) if j==0 else 1))
        mid=(l[axis]+r[axis])/2
        visit(l,r[:axis]+(mid,)+r[axis+1:],depth+1,trail+'0')
        visit(l[:axis]+(mid,)+l[axis+1:],r,depth+1,trail+'1')
    visit(*box,0,'')
    return {'index':index,'counts':dict(stats),'residual':residual,'min_weighted':str(lower) if lower is not None else None,
            'ledger_sha256':ledger.hexdigest()}


def main():
    parser=ArgumentParser();parser.add_argument('source',type=Path);parser.add_argument('output',type=Path)
    parser.add_argument('--target',default='0.837130000010');parser.add_argument('--workers',type=int,default=2)
    parser.add_argument('--samples',type=int,default=512);parser.add_argument('--depth',type=int,default=20)
    parser.add_argument('--resume',type=Path);args=parser.parse_args()
    if args.resume:
        data=json.loads(args.resume.read_text()); assert data['target']==args.target and data['frozen']==list(map(str,FROZEN))
        jobs=[(j,tuple(tuple(map(F,r)) for r in box)) for j,box in enumerate(data['residual'])]
    else:
        axes=(8,4,4);jobs=[]
        steps=[(hi-lo)/n for lo,hi,n in zip(ROOT_LO,ROOT_HI,axes)]
        for j,c in enumerate(product(*(range(n) for n in axes))):
            lo=tuple(a+q*k for a,q,k in zip(ROOT_LO,steps,c));hi=tuple(a+q for a,q in zip(lo,steps))
            jobs.append((j,(lo,hi)))
        assert sum(prod(b-a for a,b in zip(*box)) for _,box in jobs)==prod(b-a for a,b in zip(ROOT_LO,ROOT_HI))
    stats=Counter();residual=[];records=[]
    with ProcessPoolExecutor(max_workers=args.workers,initializer=initialize,initargs=(args.source,args.target,args.samples,args.depth)) as pool:
        for i,result in enumerate(pool.map(process,jobs,chunksize=1)):
            stats.update(result['counts']);residual.extend(result['residual']);records.append(result)
            if (i+1)%16==0: print(i+1,len(jobs),stats['checked'],len(residual),flush=True)
    leaves=sum(v for k,v in stats.items() if k=='taylor' or k.startswith('weights_'))+len(residual)
    assert stats['checked']==2*leaves-len(jobs)
    out={'python':platform.python_version(),'flint':flint.__version__,'numpy':np.__version__,
        'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),'checker_sha256':CHECKER_SHA,
        'target':args.target,'depth_limit':args.depth,'samples':args.samples,'frozen':list(map(str,FROZEN)),
        'root':[list(map(str,ROOT_LO)),list(map(str,ROOT_HI))],'jobs':len(jobs),'counts':dict(stats),
        'resume_sha256':sha256(args.resume.read_bytes()).hexdigest() if args.resume else None,
        'residual':residual,'records':records}
    args.output.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(args.target,dict(stats), len(residual),flush=True)


if __name__=='__main__':
    main()
