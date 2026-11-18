"""Independently integrate the n=4 packet on a shifted analytic contour."""
import argparse
import json
import os
import time
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre
from scipy.stats import qmc

from packet import A, G, W, D0, phase, required

B=3*np.sqrt(3)/8
ZSTAR=G+.5j


def integral(lam,delta,order):
    nodes,weights=roots_legendre(order)
    half_core=np.minimum(.31,8/np.sqrt(lam))
    a,b=np.full(len(lam),-3.),np.full(len(lam),3.)
    zero=np.zeros(len(lam))
    edges=np.stack((a,-G-half_core,-G+half_core,zero,G-half_core,G+half_core,b),axis=1)
    hw=(edges[:,1:]-edges[:,:-1])/2
    mid=(edges[:,1:]+edges[:,:-1])/2
    z=mid[:,:,None]+hw[:,:,None]*nodes+.5j
    kernel=np.exp(lam[:,None,None]*(phase(z,delta[:,None,None])+A-1j*B))/np.sqrt(
        1+delta[:,None,None]*z**4/2)
    return np.sqrt(lam)*np.sum(kernel*weights*hw[:,:,None],axis=(1,2))


def run(shard,repeat,exp,order):
    q=qmc.Sobol(4,scramble=True,rng=np.random.default_rng(795510+shard+15*repeat)).random_base2(exp)
    u,v,rr,zz=q.T
    d=np.exp(np.log(2.**-12)+(shard+u)/15*np.log(.75/(2.**-12)))
    lam=np.exp(np.log(16.)+v*np.log(64.))
    real_x=-3+6*rr
    horizontal_x=-3+6*zz
    dlow=required(real_x,d)
    path_excess=phase(horizontal_x+.5j,d).real+A
    cap=0.
    max_asymptotic=0.
    min_expected_modulus=np.inf
    validations=0
    within_saddle=0
    for start in range(0,len(lam),768):
        sl=slice(start,start+768)
        lamb,delta=lam[sl],d[sl]
        got=integral(lamb,delta,order)
        subset=(start+np.arange(0,len(lamb),40))
        refined=integral(lam[subset],d[subset],2*order)
        cap=max(cap,float(np.max(np.abs(refined-got[::40]))))
        assert np.all(np.isfinite(got))
        eligible=(lamb*delta>=18.)&(lamb>=128.)
        q0=1+delta[eligible]*ZSTAR**4/2
        f2=-3*ZSTAR**2-delta[eligible]/q0
        saddle=np.sqrt(2*np.pi)/(np.sqrt(q0)*np.sqrt(-f2))
        if len(saddle):
            error=np.abs(got[eligible]/saddle-1)
            max_asymptotic=max(max_asymptotic,float(np.max(error)))
            min_expected_modulus=min(min_expected_modulus,float(np.min(np.abs(got[eligible]))))
            within_saddle+=int(np.count_nonzero(error<.02))
            validations+=len(saddle)
    return dict(points=len(lam),quoted_error_required_max=float(np.max(dlow)),
                horizontal_excess_max=float(np.max(path_excess)),refined_subset=len(lam)//40,
                max_nested_quadrature_error=cap,max_eligible_relative_saddle_error=max_asymptotic,
                min_eligible_packet_modulus= float(min_expected_modulus) if validations else None,
                eligible_saddle_points=validations,eligible_within_two_percent=within_saddle)


def main():
    p=argparse.ArgumentParser()
    p.add_argument('shard',type=int)
    p.add_argument('--exponent',type=int,default=16)
    p.add_argument('--order',type=int,default=22)
    p.add_argument('--seconds',type=float,default=89.)
    p.add_argument('--max-repeats',type=int,default=160)
    p.add_argument('--out',type=Path,required=True)
    args=p.parse_args()
    start=time.monotonic()
    cpu=time.process_time()
    details=[]
    while len(details)<args.max_repeats and (not details or time.monotonic()-start<args.seconds):
        details.append(run(args.shard,len(details),args.exponent,args.order))
    out=dict(shard=args.shard,pid=os.getpid(),points=sum(r['points'] for r in details),
             elapsed=time.monotonic()-start,cpu=time.process_time()-cpu,
             seconds_target=args.seconds,order=args.order,replicas=details)
    args.out.write_text(json.dumps(out,indent=2)+"\n")


if __name__=='__main__':
    main()
