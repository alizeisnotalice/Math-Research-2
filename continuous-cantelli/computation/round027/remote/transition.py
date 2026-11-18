"""Tilted Gaussian mass and plateau widths for a quadratic minimum model."""
import argparse
import json
import math
import os
import platform
import resource
import socket
import time
from pathlib import Path

import numpy as np
import scipy
from scipy.special import log_ndtr, ndtr, ndtri_exp, roots_laguerre, roots_legendre
from scipy.stats import qmc

LOG_SQRT_2PI = .5*math.log(2*math.pi)
GL_NODES, GL_WEIGHTS = roots_legendre(64)
GL_NODES = 6*(GL_NODES+1)
GL_WEIGHTS = 6*GL_WEIGHTS
BINS = (32.,128.,512.,2048.,8192.,32768.)


def masses(a,s,b,r,lam):
    variance_inverse = 1+a*lam**2
    z = np.sqrt(variance_inverse)
    shift = lam/z
    core = np.exp(.5*shift**2)/z*(ndtr(r*z-shift)-ndtr(-r*z-shift))
    log_qm = .5*(s*lam)**2
    qb = np.exp(.5*(1-b)*lam**2)
    log_cdf_core_hi = log_ndtr(r-lam)
    correction = np.exp(log_qm+log_cdf_core_hi)*(-np.expm1(log_ndtr(-r-lam)-log_cdf_core_hi))
    center = core+qb-correction
    assert np.all((center>0)&(center<1))
    log_scale = log_qm+np.log1p(-np.exp(np.log(np.maximum(qb,1e-300))-log_qm))
    return center, log_scale


def independent_tail(q,d,s,lam,nodes,weights):
    small=q<=4
    result=np.empty_like(q)
    t=q[small,None]
    integral=np.sum(GL_WEIGHTS*np.exp(-t*GL_NODES-.5*GL_NODES**2),axis=1)
    result[small]=np.log(integral)
    t=q[~small,None]
    integral=np.sum(np.exp(-.5*(nodes/t)**2)*weights,axis=1)
    result[~small]=np.log(integral)-np.log(q[~small])
    return .5*(s-d)*(s+d)*lam**2-LOG_SQRT_2PI+result



def sample(shard,shards,exponent,replicate):
    grid=qmc.Sobol(7,scramble=True,rng=np.random.default_rng(592810+shard+shards*replicate)).random_base2(exponent)
    a,m,b,r,lam,omega,margin=grid.T
    beta=np.exp(np.log(.01)+(shard+m)/shards*np.log(2500.))
    lam=np.exp(np.log(32.)+lam*np.log(1024.))
    omega=.5+.45*(2*omega-1)*(-np.expm1(-.5*beta**2))
    return np.exp(np.log(.25)+a*np.log(32.)),beta,1.01+2*b,.025+.425*r,lam,omega,.2+.65*margin


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('shard',type=int)
    ap.add_argument('--shards',type=int,default=15)
    ap.add_argument('--exponent',type=int,default=20)
    ap.add_argument('--order',type=int,default=64)
    ap.add_argument('--replicates',type=int,default=1)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    start=time.monotonic()
    cpu_start=time.process_time()
    for replicate in range(args.replicates):
        compute(args,replicate)
    if args.replicates>1:
        receipt=dict(shard=args.shard, replicates=args.replicates,
                     pid=os.getpid(), clock_seconds=time.monotonic()-start,
                     cpu_seconds=time.process_time()-cpu_start, points=args.replicates*(2**args.exponent))
        args.out.write_text(json.dumps(receipt,indent=2)+"\n")


def compute(args,replicate):
    start=time.monotonic()
    nodes,weights=roots_laguerre(args.order)
    a,beta,b,r,lam,omega,margin=sample(args.shard,args.shards,args.exponent,replicate)
    count=len(a)
    s=beta/lam
    c0,log_scale=masses(a,s,b,r,lam)
    log_p_star=np.log1p(-c0)-log_scale
    d_sym=-ndtri_exp(log_p_star-math.log(2.))/lam
    d1=-ndtri_exp(log_p_star+np.log(omega))/lam
    d2=-ndtri_exp(log_p_star+np.log1p(-omega))/lam
    assert np.all(d1+d2>=2*d_sym-1e-14)
    assert np.all((1-d1)*lam>r)
    assert np.all((d1>0)&(d2>0))
    log_q0_arg=-.5*beta**2-math.log(2.)
    q0=-ndtri_exp(log_q0_arg)
    c_coeff=np.exp(.5/a)/np.sqrt(a)
    expected_shift=c_coeff*np.exp(log_q0_arg+.5*q0*q0+LOG_SQRT_2PI)
    d_sym_asym=(q0+expected_shift/lam)/lam
    correction=(d_sym*lam-q0)*lam-expected_shift
    # Narrowest symmetric finite witness for budget margin in the model.
    witness_log_mass=np.log(margin)+log_p_star
    dw=-ndtri_exp(witness_log_mass-math.log(2.))/lam
    assert np.all((1-dw)*lam>r)
    target=c0+margin*(1-c0)
    max_rel=0.
    max_cross=0.
    maximum_location={}
    final_max=0.
    for lo in range(0,count,8192):
        i=slice(lo,lo+8192)
        verified=c0[i]+np.exp(log_scale[i]+np.logaddexp(log_ndtr(-lam[i]*d1[i]),log_ndtr(-lam[i]*d2[i])))
        assert np.all(np.abs(verified-1.)<5e-6)
        f_test=c0[i]+2*np.exp(log_scale[i]+log_ndtr(-lam[i]*dw[i]))
        assert np.all(np.abs(f_test-target[i])<5e-6)
        logs1=independent_tail(lam[i]*d1[i],d1[i],s[i],lam[i],nodes,weights)
        logs2=independent_tail(lam[i]*d2[i],d2[i],s[i],lam[i],nodes,weights)
        analytical1=log_ndtr(-lam[i]*d1[i]) + .5*beta[i]**2
        analytical2=log_ndtr(-lam[i]*d2[i]) + .5*beta[i]**2
        array_err=np.maximum(np.abs(logs1-analytical1),np.abs(logs2-analytical2))
        if float(max(array_err)) > max_cross:
            j=lo+int(np.argmax(array_err))
            max_cross=float(max(array_err))
            maximum_location=dict(a=float(a[j]),beta=float(beta[j]),b=float(b[j]),r=float(r[j]),
                                  lam=float(lam[j]),omega=float(omega[j]))
        max_rel=max(max_rel,float(max(np.abs(verified-1.))))
        final_max=max(final_max,float(max(np.abs(f_test-target[i]))))
    bins=[]
    for low,high in zip(BINS[:-1],BINS[1:]):
        idx=(lam>=low)&(lam<high)
        bins.append(dict(lam=[low,high],samples=int(np.count_nonzero(idx)),
                         max_abs_width_error=float(max(np.abs(d_sym[idx]-d_sym_asym[idx]))),
                         max_abs_action_shift_error=float(max(np.abs(correction[idx]))),
                         max_width_correction=float(max((d_sym*lam-q0)[idx]))))
    out=dict(shard=args.shard, replicate=replicate, host=socket.gethostname(), machine=platform.machine(), python=platform.python_version(),
             numpy=np.__version__,scipy=scipy.__version__, order=args.order, points=count, pid=os.getpid(),
             clock_seconds=time.monotonic()-start,cpu_seconds=time.process_time(),
             max_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,beta_range=[float(min(beta)),float(max(beta))],
             max_float64_threshold_residual=max_rel, max_finite_witness_residual=final_max,
             max_independent_log_action_error=max_cross, worst_independent_case=maximum_location, by_lambda=bins)
    result = args.out if args.replicates==1 else args.out.with_name(f"{args.out.stem}-r{replicate}.json")
    result.write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps(dict(shard=args.shard,points=count,secs=out['clock_seconds'],pid=os.getpid())))


if __name__=='__main__':
    main()
