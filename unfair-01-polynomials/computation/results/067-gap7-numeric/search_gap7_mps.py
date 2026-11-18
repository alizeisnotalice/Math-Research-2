"""Targeted numerical search for a finite extension of the exact C7-gap prefix."""
from argparse import ArgumentParser
from fractions import Fraction
import gc
import hashlib
import json
import os
from pathlib import Path
import sys
import time

os.environ['PYTORCH_ENABLE_MPS_FALLBACK']='0'
import numpy as np
import sympy as sp
import torch


def atomic(path, data):
    tmp=path.with_suffix(path.suffix+'.partial')
    tmp.write_text(json.dumps(data,indent=2)+'\n')
    tmp.replace(path)


def prefix():
    t,x=sp.symbols('t x')
    H=sp.Poly(t**6-4*t**5+6*t**4-7*t**3+4*t**2-2*t+1,t)
    assert H.count_roots(0,1)==1 and sp.Poly(H,modulus=2).is_irreducible
    m=1-t+t*t
    n=1-t*m
    p=1-(1-t)*n
    v=t*t*m
    L=[1,t,0,n,p,0,0,0]
    R=[1,1-t,m,0,0,v,0,0]
    conv=sp.Poly(sp.expand(sum(z*x**i for i,z in enumerate(L))*sum(z*x**i for i,z in enumerate(R))),x)
    assert [sp.rem(conv.nth(i),H.as_expr(),t) for i in range(8)]==[1]*7+[0]
    assert sp.expand(conv.nth(8)-n*v)==0
    lo,hi=Fraction(7,10),Fraction(3,4)
    def evaluate(z):
        out=Fraction(0)
        for coeff in [1,-4,6,-7,4,-2,1]:out=out*z+coeff
        return out
    assert evaluate(lo)>0>evaluate(hi)
    for _ in range(66):
        mid=(lo+hi)/2
        if evaluate(mid)>0:
            lo=mid
        else:
            hi=mid
    tf=float((lo+hi)/2)
    values=lambda expr:float(sp.sympify(expr).subs(t,tf))
    return [values(x) for x in L],[values(x) for x in R],values(1-n*v),{
        'minpoly_ascending':[1,-2,4,-7,6,-4,1],
        'isolation_interval':[str(lo),str(hi)],
        'left':list(map(str,L)),'right':list(map(str,R)),
        'degree8_remaining_budget':str(1-n*v),
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }


def main():
    ap=ArgumentParser()
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--m',type=int,default=16)
    ap.add_argument('--n',type=int,default=24)
    ap.add_argument('--starts',type=int,default=4096)
    ap.add_argument('--steps',type=int,default=16000)
    ap.add_argument('--restart',type=int,default=2000)
    ap.add_argument('--seconds',type=float,default=14400)
    ap.add_argument('--seed',type=int,default=2026111807)
    a=ap.parse_args()
    assert __debug__ and Path('/Volumes/data').is_mount()
    assert min(a.m,a.n)>=10 and a.starts>=16 and a.steps>0 and a.restart>0
    a.out.mkdir(parents=True)
    assert torch.backends.mps.is_available()
    torch.set_num_threads(1)
    torch.manual_seed(a.seed)
    left,right,budget,prefix_data=prefix()
    assert 0<budget<1
    params={k:str(v) if isinstance(v,Path) else v for k,v in vars(a).items()}
    atomic(a.out/'input.json',{'parameters':params,'prefix':prefix_data,'python':sys.version,
        'torch':torch.__version__,'numpy':np.__version__,'sympy':sp.__version__,'device':'mps','dtype':'float32',
        'scope':'Numerical optimization for an exact C7-gap starting prefix. Samples alone have no exclusion credit.'})
    device=torch.device('mps');k=a.starts;m,n=a.m,a.n
    la=torch.tensor(left,dtype=torch.float32,device=device).expand(k,-1)
    lb=torch.tensor(right,dtype=torch.float32,device=device).expand(k,-1)
    terminal=torch.ones((k,1),device=device)
    incidence=np.zeros(((m+1)*(n+1),m+n+1),dtype=np.float32)
    for i in range(m+1):
        for j in range(n+1):incidence[i*(n+1)+j,i+j]=1
    select=torch.tensor(incidence,device=device)
    del incidence
    lock_a=sorted({j-d for j in (m,n) for d in (1,2,5) if 0<=j-d<=m})
    lock_b=sorted({j-d for j in (m,n) for d in (1,3,4) if 0<=j-d<=n})
    assert all(12<=i<m for i in lock_a) and all(12<=i<n for i in lock_b)
    keep_a=torch.tensor([int(i not in lock_a) for i in range(m+1)],device=device)
    keep_b=torch.tensor([int(i not in lock_b) for i in range(n+1)],device=device)
    atomic(a.out/'saturation.json',{'A_zero_positions':lock_a,'B_zero_positions':lock_b,
        'fixed_unit_product_degrees':list(range(12,31)),
        'score':'0.7*mean_squared_binary_residue + 0.3*maximum_squared_binary_residue'})
    theta=torch.nn.Parameter(torch.zeros((k,4),device=device))
    free=torch.nn.Parameter(torch.zeros((k,m+n-24),device=device))
    opt=torch.optim.Adam([theta,free],lr=.027)
    history=[]
    windows=[]
    best=[]
    start=time.monotonic()
    completed=0
    status='failed'
    try:
        for step in range(1,a.steps+1):
            if time.monotonic()-start>=a.seconds:
                status='partial_time';break
            if (step-1)%a.restart==0:
                if best:windows.append({'last_step':completed,'cases':best})
                with torch.no_grad():
                    theta.uniform_(0,1)
                    theta[:k//8].zero_()
                    theta[k//8:k//4].fill_(1)
                    free.uniform_(0,.22)
                opt.state.clear()
            opt.zero_grad(set_to_none=True)
            a8,b8=budget*theta[:,0:1],budget*(1-theta[:,0:1])
            a1,b1=la[:,1:2],lb[:,1:2]
            b2,a3,a4,b5=lb[:,2:3],la[:,3:4],la[:,4:5],lb[:,5:6]
            budget9=1-a4*b5-a1*b8-b1*a8
            a9,b9=budget9*theta[:,1:2],budget9*(1-theta[:,1:2])
            budget10=1-a1*b9-b1*a9-b2*a8
            a10,b10=budget10*theta[:,2:3],budget10*(1-theta[:,2:3])
            budget11=(1-a1*b10-b1*a10-b2*a9-a3*b8).clamp(min=0)
            a11,b11=budget11*theta[:,3:4],budget11*(1-theta[:,3:4])
            aa=torch.cat((la,a8,a9,a10,a11,free[:,:m-12],terminal),dim=1)
            bb=torch.cat((lb,b8,b9,b10,b11,free[:,m-12:],terminal),dim=1)
            aa=aa*keep_a
            bb=bb*keep_b
            cv=(aa[:,:,None]*bb[:,None,:]).reshape(k,-1)@select
            assert cv.device.type=='mps'
            target=(cv[:,9:].detach()>=.5).to(cv.dtype)
            target[:,3:22]=1
            residue=cv[:,9:]-target
            squared=residue.square()
            scores=0.7*squared.mean(dim=1)+0.3*squared.amax(dim=1)
            scores.sum().backward()
            opt.step()
            with torch.no_grad():
                theta.clamp_(0,1);free.clamp_(0,1)
            completed=step
            if step==1 or step%200==0 or step==a.steps:
                with torch.no_grad():
                    ix=torch.topk(scores.detach(),16,largest=False).indices
                    lefts=aa.detach()[ix].cpu().numpy();rights=bb.detach()[ix].cpu().numpy()
                    best=[]
                    for av,bv in zip(lefts,rights,strict=True):
                        conv=np.convolve(av.astype(np.float64),bv.astype(np.float64))
                        rt=(conv>=.5).astype(np.float64)
                        assert list(rt[:9])==[1.]*7+[0.,1.]
                        best.append({'A':av.tolist(),'B':bv.tolist(),'residual_linf':float(np.max(np.abs(conv-rt)))})
                    best.sort(key=lambda z:z['residual_linf'])
                    record={'step':step,'lowest_residual':best[0]['residual_linf'],'device':str(cv.device),
                            'driver_bytes':torch.mps.driver_allocated_memory(),'seconds':time.monotonic()-start}
                    history.append(record)
                    atomic(a.out/'checkpoint.json',{'history':history,'current':best,'completed_windows':windows})
                    print(json.dumps(record),flush=True)
        else:
            status='finished_numeric_search'
        atomic(a.out/'result.json',{'status':status,'steps':completed,'history':history,'completed_windows':windows,'current':best})
    finally:
        del opt,theta,free,la,lb,select,terminal,keep_a,keep_b
        if completed:
            del aa,bb,cv,residue,scores,squared
        gc.collect();torch.mps.synchronize();torch.mps.empty_cache()
        atomic(a.out/'release.json',{'cleanup_called':True,'terminal_status':status,'step':completed,'pid':os.getpid()})


if __name__=='__main__':
    main()
