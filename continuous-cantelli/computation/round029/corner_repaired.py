"""Minimum-width Gaussian plateau constrained by the origin-core boundary."""
import argparse
import json
import math
import os
import platform
import time
from pathlib import Path

import numpy as np
import scipy
from scipy.special import log_ndtr, ndtr, ndtri_exp, roots_laguerre
from scipy.stats import qmc

LOG2=math.log(2)
LOGSQRT2PI=.5*math.log(2*math.pi)


def right_quantile(log_tail):
    q=-ndtri_exp(log_tail)
    # One asymptotically stable log-CDF Newton update removes the
    # near-core branch misplacement caused by the raw inverse residual.
    value=log_ndtr(-q)
    hazard=np.exp(-.5*q*q-LOGSQRT2PI-value)
    return q+(value-log_tail)/hazard


def run(shard,repeat,count,order):
    x=qmc.Sobol(6,scramble=True,rng=np.random.default_rng(209651+shard+15*repeat)).random_base2(count)
    aa,bb,rr,ll,pp,tt=x.T
    a=.5+7.5*aa
    b=1.01+2*bb
    r=.30+.95*rr
    lam=32*np.exp(np.log(128.)*ll)
    phase=2*r*np.exp2(-2+4*(shard+pp)/15)
    tau=-2+4*tt
    near=2*r-2*(np.log(lam)+tau)/lam
    beta=np.where(np.arange(len(x))%2==0,phase,near)
    assert np.all((beta>0)&(beta<lam))
    z=np.sqrt(1+a*lam**2)
    j=np.exp(.5*(lam/z)**2)/z*(ndtr(r*z-lam/z)-ndtr(-r*z-lam/z))
    lm=.5*lam*(lam-beta)
    lb=.5*(1-b)*lam**2
    qb=np.exp(lb)
    lq=lm+np.log1p(-np.exp(lb-lm))
    far_left=np.exp(lm+log_ndtr(-lam-r))
    ltouch=log_ndtr(r-lam)
    ccon=j+qb*(-np.expm1(ltouch))+far_left
    l_remaining=np.log1p(-ccon)
    constrained=lq+ltouch>=l_remaining
    critical=right_quantile(l_remaining-lq)
    sym=~constrained
    core_scale=np.exp(lm[sym]+ltouch[sym])*(-np.expm1(log_ndtr(-lam[sym]-r[sym])-ltouch[sym]))
    baseline=j[sym]+qb[sym]-core_scale
    critical[sym]=right_quantile(np.log1p(-baseline)-lq[sym]-LOG2)
    left=np.where(constrained,r,lam-critical)
    assert np.all(left>=r-1e-10), float(np.min(left-r))
    left=np.maximum(left,r)
    right=lam+critical
    assert np.all(right>lam)
    lcdf_a=log_ndtr(left-lam)
    inner_left=np.zeros_like(left)
    inner_left[sym]=np.exp(lm[sym]+lcdf_a[sym])*(-np.expm1(np.minimum(0,ltouch[sym]-lcdf_a[sym])))
    rtail=np.exp(lm+log_ndtr(-critical))
    f=j+far_left+inner_left+rtail+qb*(1-np.exp(lcdf_a)-np.exp(log_ndtr(-critical)))
    assert np.all(np.abs(f-1)<2e-8), float(np.max(np.abs(f-1)))
    nodes,weights=roots_laguerre(order)
    max_independent=0.
    for start in range(0,len(x),8192):
        sl=slice(start,start+8192)
        q=critical[sl]
        quad_raw=np.sum(np.exp(-.5*(nodes/q[:,None])**2)*weights,axis=1)
        independent=lm[sl]-.5*q*q-np.log(q)-LOGSQRT2PI+np.log(quad_raw)
        expected=lm[sl]+log_ndtr(-q)
        max_independent=max(max_independent,float(np.max(np.abs(independent-expected))))
    # The mass at the touching threshold tends to
    # exp(tau-r²/2)/sqrt(2*pi) in the log-narrow transition family.
    is_near=np.arange(len(x))%2==1
    limits=tau-r*r/2-LOGSQRT2PI
    observed=lq+ltouch
    return dict(shard=shard,repeat=repeat,point_count=len(x),constrained=int(np.count_nonzero(constrained)),
                near_count=int(np.count_nonzero(is_near)),
                phase_counts=[int(np.count_nonzero(constrained & ~is_near)),int(np.count_nonzero(constrained & is_near))],
                max_threshold_residual=float(np.max(np.abs(f-1))),max_quadrature_log_difference=max_independent,
                max_phase_error_high_lam=float(np.max(np.abs(observed[is_near & (lam>=1024)]-
                                                            limits[is_near & (lam>=1024)]))))


def main():
    p=argparse.ArgumentParser()
    p.add_argument('shard',type=int)
    p.add_argument('--exponent',type=int,default=20)
    p.add_argument('--replicates',type=int,default=16)
    p.add_argument('--order',type=int,default=112)
    p.add_argument('--out',type=Path,required=True)
    args=p.parse_args()
    elapsed=time.monotonic()
    cpu=time.process_time()
    results=[]
    for i in range(args.replicates):
        results.append(run(args.shard,i,args.exponent,args.order))
    out=dict(shard=args.shard,pid=os.getpid(),points=sum(r['point_count'] for r in results),
             worker_wall_seconds=time.monotonic()-elapsed,worker_cpu_seconds=time.process_time()-cpu,
             python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,results=results)
    args.out.write_text(json.dumps(out,indent=2)+"\n")


if __name__=='__main__':
    main()
