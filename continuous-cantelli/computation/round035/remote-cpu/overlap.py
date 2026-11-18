"""Gaussian dual overlap of a critical-origin core and a disjoint tail patch."""
import argparse
import json
import math
import os
import platform
import time
from pathlib import Path

import numpy as np
import scipy
from scipy.special import ndtr, roots_legendre
from scipy.stats import qmc

SQRT2PI=math.sqrt(2*math.pi)
LOGSQRT2=.5*math.log(2)
N8,W8=roots_legendre(8)
N20,W20=roots_legendre(20)


def phi(s):
    return np.exp(-.5*s*s)/SQRT2PI


def core(L,t,j,q,obs,half,order,direct):
    n,w=roots_legendre(order)
    y=L[:,None]*(3.25+2.75*n)
    ell=L[:,None]-np.log(y)
    d=ell+.5*np.log(ell)+LOGSQRT2+(j[:,None]*np.log(ell)+q[:,None])/ell
    assert np.all(d>0)
    scale=np.sqrt(2*d)
    z=((1+obs[:,None]/t[:,None])*y-2*d)/scale
    span=half[:,None]*y/(t[:,None]*scale)
    if direct:
        averaged=np.sum(W8*np.exp(-.5*(z[:,:,None]+span[:,:,None]*N8)**2),axis=2)/2
    else:
        averaged=np.empty(z.shape)
        narrow=span<.005
        zn,ss=z[narrow],span[narrow]**2
        averaged[narrow]=np.exp(-.5*zn*zn)*(1+(zn*zn-1)*ss/6+(zn**4-6*zn*zn+3)*ss**2/120)
        zz,hh=z[~narrow],span[~narrow]
        averaged[~narrow]=(ndtr(zz+hh)-ndtr(zz-hh))*SQRT2PI/(2*hh)
    return 2.75*L*np.sum(w*averaged*np.exp(d-L[:,None])*phi(y/t[:,None]),axis=1)


def tail(t,k,offset,half_patch,obs,half_obs,direct):
    s=obs[:,None]+half_obs[:,None]*N20
    amp=np.exp(-.5*k[:,None]*(1+s/t[:,None])**2)
    if direct:
        z=offset[:,None]+half_patch[:,None]*N20
        prob=half_patch[:,None]*np.sum(W20*phi(z[:,None,:]-s[:,:,None]),axis=2)
    else:
        prob=ndtr((offset+half_patch)[:,None]-s)-ndtr((offset-half_patch)[:,None]-s)
    return np.sum(W20*amp*prob,axis=1)/2


def evaluate(pars,order,direct=False):
    L,j,q,k,offset,width,center,half=pars.T
    t=np.exp(L)
    # Source core y in [L/2,6L] and the remote source patch are disjoint.
    assert np.all(6*L/t < t+offset-width)
    c=core(L,t,j,q,center,half,order,direct)
    a=tail(t,k,offset,width,center,half,direct)
    p0=core(L,t,j,q,np.zeros_like(t),np.zeros_like(t),order,True)
    t0=tail(t,k,offset,width,np.zeros_like(t),np.zeros_like(t),False)
    return c,a,p0+t0


def run(shard,repeat,exponent,order):
    grid=qmc.Sobol(8,scramble=True,rng=np.random.default_rng(811704+shard+15*repeat)).random_base2(exponent)
    l,j,q,k,u,w,c,h=grid.T
    pars=np.column_stack((7+17*(shard+l)/15,-.75+1.4*j,-2.5+4*q,-.5+7*k,
                          -6+12*u,.05+2.4*w,-6+12*c,.04+3.3*h))
    best=[]
    count_upper=0
    count_dual_improvement=0
    max_dis=0.
    for start in range(0,len(grid),768):
        sample=pars[start:start+768]
        pc,pt,point=evaluate(sample,order)
        score=pc+pt
        count_upper+=int(np.count_nonzero(score>1))
        qualified=(score>1)&(point<=1)
        count_dual_improvement+=int(np.count_nonzero(qualified))
        arg=np.argpartition(score,-min(2,len(score)))[-min(2,len(score)):]
        qarg=np.flatnonzero(qualified)
        if len(qarg):
            arg=np.union1d(arg,np.array([qarg[np.argmax(score[qarg])]]))
        compare=np.union1d(arg,np.arange(0,len(sample),64))
        dc,dt,dpoint=evaluate(sample[compare],2*order,True)
        max_dis=max(max_dis,float(np.max(np.abs(dc+dt-score[compare]))))
        for i,row_c,row_t,reference in zip(compare,dc,dt,dpoint,strict=True):
            if i in arg:
                best.append(dict(parameters=sample[i].tolist(),core=float(pc[i]),tail=float(pt[i]),
                                 dual=float(score[i]),point_at_t=float(point[i]),
                                 independent_dual=float(row_c+row_t),independent_point_at_t=float(reference)))
        if len(best)>50:
            strong=sorted(best,key=lambda r:r['dual'],reverse=True)[:6]
            improvement=sorted([r for r in best if r['point_at_t']<=1],key=lambda r:r['dual'],reverse=True)[:6]
            best=strong+improvement
    return dict(shard=shard,repeat=repeat,points=len(grid),positive=count_upper,
                positive_despite_nominal_point=count_dual_improvement,max_crosscheck=max_dis,top=best)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('shard',type=int)
    ap.add_argument('--exponent',type=int,default=19)
    ap.add_argument('--order',type=int,default=64)
    ap.add_argument('--replicates',type=int,default=3)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    wall,cpu=time.monotonic(),time.process_time()
    results=[]
    for i in range(args.replicates):
        results.append(run(args.shard,i,args.exponent,args.order))
    out=dict(shard=args.shard,pid=os.getpid(),points=sum(r['points'] for r in results),
             wall_seconds=time.monotonic()-wall,cpu_seconds=time.process_time()-cpu,
             python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
             parameter_order=['L','j','q','k','x_patch_mid_offset','x_patch_halfwidth',
                              'dual_window_mid_offset','dual_window_halfwidth'],results=results)
    args.out.write_text(json.dumps(out,indent=2)+"\n")


if __name__=='__main__':
    main()
