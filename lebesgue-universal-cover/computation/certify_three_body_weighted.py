"""All-angle rational convex weights for overlapping moving inner fans."""
from argparse import ArgumentParser
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec,spec_from_file_location
from itertools import product
import json
from math import floor
from pathlib import Path

from flint import arb
import numpy as np
from scipy.optimize import linprog, OptimizeWarning
import warnings

from certify_three_body_refinement import CHECKER_SHA,SEED_SHA,ROOT_LO,ROOT_HI,split,determinant

SCALE=10**12
SIGNS=tuple(product((-1,1),repeat=4))


def area_coefficients(runs,cs,sn):
    """Twice area = c0 + linear dot (t3x,t3y,t5x,t5y) + c5 det(t3,t5)."""
    def rot(body,q):
        return q if body!=2 else (cs*q[0]-sn*q[1],sn*q[0]+cs*q[1])
    out=[arb(0) for _ in range(6)]
    def shift(body,x,y):
        if body:
            pos=2*body-1;out[pos]+=x;out[pos+1]+=y
    for k,(i,s,e,intrinsic) in enumerate(runs):
        start,end=rot(i,s),rot(i,e)
        j,ns=runs[(k+1)%len(runs)][:2];v=rot(j,ns)
        out[0]+=intrinsic+determinant(end,v)
        shift(i,end[1]-start[1]+v[1],start[0]-end[0]-v[0]);shift(j,-end[1],end[0])
        if (i,j)==(1,2):out[5]+=1
        if (i,j)==(2,1):out[5]-=1
    return out


def corner_row(runs,pa,h):
    f0=area_coefficients(runs,1,0); fpi=area_coefficients(runs,-1,0);f90=area_coefficients(runs,0,1)
    a=[(v+w)/2 for v,w in zip(f0,fpi)]; b=[v-w for v,w in zip(f0,a)]; c=[v-w for v,w in zip(f90,a)]
    mid,rad=pa[2][:2];cs=mid.cos();sn=mid.sin();
    fc=[v+w*cs+x*sn for v,w,x in zip(a,b,c)];fp=[-w*sn+x*cs for w,x in zip(b,c)]
    def at(coeff,t):return sum((v*w for v,w in zip(coeff,t)),arb(0))
    center=(*pa[1][2],*pa[2][2]);width=tuple(arb(float(h[k])) for k in (0,1,3,4));out=[]
    for signs in SIGNS:
        x,y,u,v=[z+s*w for z,s,w in zip(center,signs,width)];basis=(1,x,y,u,v,x*v-y*u)
        # A subinterval lower bound from midpoint derivative and a full-interval curvature bound.
        lb=(at(fc,basis)-abs(at(fp,basis))*rad-(abs(at(b,basis))+abs(at(c,basis)))*rad**2/2)/2
        value=floor(float(lb.lower())*SCALE)-1
        assert lb>arb(value)/SCALE
        out.append(value)
    return out


def dual(rows,target):
    """Discover weights in floating point; accept only strict integer corner inequalities."""
    if any(min(row)>target for row in rows):
        j=next(j for j,row in enumerate(rows) if min(row)>target);return [(j,1)]
    if len(rows)<2:return None
    if any(max(row[v] for row in rows)<=target for v in range(16)):return None
    # At most one new row was appended since the last dual search.
    b=rows[-1]
    for j,a in enumerate(rows[:-1]):
        lo=Q(0);hi=Q(1)
        for x,y in zip(a,b):
            if x==y:
                if y<=target:break
            elif x>y:lo=max(lo,Q(target-y,x-y))
            else:hi=min(hi,Q(target-y,x-y))
            if lo>=hi:break
        else:
            w=(lo+hi)/2;weights=[(j,w.numerator),(len(rows)-1,w.denominator-w.numerator)]
            if valid_weights(rows,weights,target):return weights
    if len(rows)<3:return None
    # Scale margins before the heuristic LP; rational row checks absorb solver tolerances.
    gaps=(np.array(rows,dtype=np.int64)-target)/1e6
    m=len(rows);mat=np.column_stack([-gaps.T,np.ones(16)])
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=OptimizeWarning, message='Unrecognized options detected.*')
        res=linprog(np.r_[np.zeros(m),-1.], A_ub=mat,b_ub=np.zeros(16),
                    A_eq=[np.r_[np.ones(m),0.]],b_eq=[1.],bounds=[(0,None)]*m+[(None,None)],
                    method='highs',options={'threads':1})
    if not res.success:return None
    weights=[(j,round(float(v)*10**8)) for j,v in enumerate(res.x[:-1]) if v>0.]
    return weights if valid_weights(rows,weights,target) else None


def valid_weights(rows,weights,target):
    total=sum(a for _,a in weights)
    return total>0 and all(a>=0 for _,a in weights) and all(sum(a*rows[j][v] for j,a in weights)>target*total for v in range(16))


def init_worker(source,target,max_depth):
    global checker, threshold_int, max_extra
    file=Path(source)/'checker/lbcheck.py';assert sha256(file.read_bytes()).hexdigest()==CHECKER_SHA
    spec=spec_from_file_location('tomar',file);lib=module_from_spec(spec);spec.loader.exec_module(lib)
    checker=lib.Checker([3,5],float(target),384,3e-6);checker.tau=arb(target)
    ti=Q(target)*SCALE;assert ti.denominator==1;threshold_int=ti.numerator;max_extra=max_depth


def solve(job):
    part,idx,box=job;stats=Counter();unresolved=[];log=sha256();smallest=None;deepest=0
    def visit(l,r,dep,path):
        nonlocal smallest,deepest
        stats['checked']+=1;deepest=max(deepest,dep)
        cp=[(a+b)/2 for a,b in zip(l,r)];wp=[(b-a)/2 for a,b in zip(l,r)];assert all(Q(float(z))==z for z in cp+wp)
        c=np.array(cp,dtype=float);h=np.array(wp,dtype=float);rows=[];shifts=[];seen=set()
        def attempt(shift):
            nonlocal smallest
            shift=tuple(shift)
            if shift in seen:return False
            seen.add(shift);runs,pa=checker.verified_runs(c,h,c+h*np.array(shift))
            if runs is None:return False
            row=corner_row(runs,pa,h);rows.append(row);shifts.append(shift)
            weights=dual(rows,threshold_int)
            if weights is None:return False
            assert valid_weights(rows,weights,threshold_int)
            wpairs=[(j,a) for j,a in weights if a];den=sum(a for _,a in wpairs)
            lower=min(Q(sum(a*rows[j][v] for j,a in wpairs),den*SCALE) for v in range(16));assert lower>Q(threshold_int,SCALE)
            if smallest is None or lower<smallest:smallest=lower
            stats['weights_'+str(len(wpairs))]+=1
            log.update(json.dumps([path,[(shifts[j],a) for j,a in wpairs],str(lower)],separators=(',',':')).encode());return True
        if attempt((0,0,0,0,0)):return
        index=min(range(16),key=lambda v:rows[0][v]) if rows else 0
        x,y,u,v=SIGNS[index]
        proposals=[(x,y,0,u,v),(0,0,-1,0,0),(0,0,1,0,0),(x,y,-1,u,v),(x,y,1,u,v),(-x,-y,0,-u,-v),(x,y,0,0,0),(0,0,0,u,v)]
        for shift in proposals:
            if attempt(shift):return
        if dep==max_extra:
            unresolved.append([[str(z) for z in l],[str(z) for z in r]]);return
        axis=max(range(5),key=lambda j:wp[j]*(Q(22,100) if j==2 else 1))
        for k,(ll,rr) in enumerate(split(l,r,axis)):visit(ll,rr,dep+1,path+str(k))
    visit(*box,0,'')
    return {'part':part,'index':idx,'counts':dict(stats),'unresolved':unresolved,'min_bound':str(smallest) if smallest is not None else None,'depth':deepest,'terminal_evidence_sha256':log.hexdigest()}


def main():
    ap=ArgumentParser();ap.add_argument('source');ap.add_argument('stage',type=Path);ap.add_argument('--workers',type=int,default=4);ap.add_argument('--max-depth',type=int,default=15);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    raw=a.stage.read_bytes();s=json.loads(raw);assert s['checker_sha256']==CHECKER_SHA and s['seed_sha256']==SEED_SHA and s['root']==[[str(x) for x in row] for row in [ROOT_LO,ROOT_HI]]
    assert [p['part'] for p in s['parts']]==list(range(16))
    jobs=[(p['part'],j,(tuple(map(Q,l)),tuple(map(Q,r)))) for p in s['parts'] for j,(l,r) in enumerate(p['unresolved'])];assert len(jobs)==s['unresolved']
    counts=Counter();parts=[];bad=0
    with ProcessPoolExecutor(max_workers=a.workers,initializer=init_worker,initargs=(a.source,s['target'],a.max_depth)) as ex:
        for k,p in enumerate(ex.map(solve,jobs,chunksize=1)):
            counts.update(p['counts']);parts.append(p);bad+=len(p['unresolved'])
            if (k+1)%200==0:print(f'{k+1}/{len(jobs)} initial residual cells; {counts["checked"]} visits; {bad} residual',flush=True)
    leaves=sum(v for k,v in counts.items() if k.startswith('weights_'))+bad
    assert counts['checked']==2*leaves-len(jobs)
    bounds=[Q(p['min_bound']) for p in parts if p['min_bound'] is not None]
    out={'target':s['target'],'source_stage_sha256':sha256(raw).hexdigest(),'completed':bad==0,'residual_input':len(jobs),'residual_output':bad,'counts':dict(counts),'max_depth':a.max_depth,'min_weighted_bound':str(min(bounds)),'parts':parts}
    a.output.write_text(json.dumps(out,separators=(',',':'))+'\n');print({k:v for k,v in out.items() if k!='parts'},flush=True)

if __name__=='__main__':main()
