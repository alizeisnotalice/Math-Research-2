"""Real folds of an explicit smooth divergent quartic jet."""
import argparse
import json
import os
import platform
import time
from pathlib import Path

import numpy as np
import scipy
from scipy.special import roots_laguerre, roots_legendre
from scipy.stats import qmc

SN, SW = roots_legendre(128)
SN = (SN + 1)/2
SW = SW/2
ALT_S = SN/(1-SN)
ALT_W = SW*np.exp(-ALT_S)/(1-SN)**2


def solver(r,nodes,weights):
    x=.5+r/4
    for _ in range(4):
        v=1+(r*x*x)[:,None]*nodes
        inv=1/v
        a=np.sum(inv**2*weights,axis=1)
        ap=np.sum(inv**3*nodes*weights,axis=1)
        x -= (2*x*a-1)/(2*a-8*r*x*x*ap)
    v=1+(r*x*x)[:,None]*nodes
    inv=1/v
    h=np.sum(inv*weights,axis=1)
    a=np.sum(inv**2*weights,axis=1)
    return x, -np.log(4*(x-x*x*h)), abs(2*x*a-1)


def scan(sid,rep,exp,nodes,weights):
    draw=qmc.Sobol(4,scramble=True,rng=np.random.default_rng(712030+sid+15*rep)).random_base2(exp)
    u,v,c,d=draw.T
    gamma=np.exp2(-3+6*(sid+u)/15)
    b=np.exp2(3+5*v)
    n=np.exp2(5+15*c)
    var=np.exp2(-2+4*d)
    r=gamma/b**2
    maxima=np.zeros(5)
    refined_count=0
    roots_min=1
    roots_max=0
    for lo in range(0,len(r),1536):
        i=slice(lo,lo+1536)
        x,log_ratio,residual=solver(r[i],nodes,weights)
        assert np.all((.5<=x)&(x<=.6))
        target_log=-r[i]/4
        approximate_error=np.abs((log_ratio-target_log)/r[i]**2)
        analytic_error=np.abs(n[i]*log_ratio+gamma[i]/(4*var[i]))
        shell=np.abs(b[i]/np.sqrt(var[i]*n[i])-1)<.04
        ids=np.arange(0,len(x),29)
        xa,la,ra=solver(r[i][ids],ALT_S,ALT_W)
        local=[float(max(residual)),float(max(abs(x[ids]-xa))),
               float(max(abs(log_ratio[ids]-la))),float(max(approximate_error)),
               float(max(analytic_error[shell])) if np.any(shell) else 0.]
        maxima=np.maximum(maxima,local)
        refined_count+=len(ids)
        roots_min=min(roots_min,float(min(x)))
        roots_max=max(roots_max,float(max(x)))
    assert maxima[0]<1e-13 and maxima[1]<1e-11 and maxima[2]<1e-11
    return dict(points=len(r),refined=refined_count,max_stationarity=float(maxima[0]),
                max_independent_root_difference=float(maxima[1]),max_independent_log_ratio_difference=float(maxima[2]),
                max_second_order_scaled_remainder=float(maxima[3]),
                max_shell_error_4percent=float(maxima[4]),range_scaled_folds=[roots_min,roots_max])


def main():
    p=argparse.ArgumentParser()
    p.add_argument('shard',type=int)
    p.add_argument('--exponent',type=int,default=14)
    p.add_argument('--seconds',type=float,default=46.)
    p.add_argument('--out',type=Path,required=True)
    args=p.parse_args()
    st,pc=time.monotonic(),time.process_time()
    nodes,weights=roots_laguerre(48)
    rows=[]
    while len(rows)<4096 and (not rows or time.monotonic()-st<args.seconds):
        rows.append(scan(args.shard,len(rows),args.exponent,nodes,weights))
    out=dict(shard=args.shard,pid=os.getpid(),points=sum(x['points'] for x in rows),
             elapsed=time.monotonic()-st,cpu=time.process_time()-pc,python=platform.python_version(),
             numpy=np.__version__,scipy=scipy.__version__,results=rows)
    args.out.write_text(json.dumps(out,indent=2)+"\n")


if __name__=='__main__':
    main()
