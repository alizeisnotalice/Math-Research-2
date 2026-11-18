"""Scan explicit n=4 Gaussian packet action and a short real input contour."""
import argparse
import json
import math
import os
import platform
import time
from pathlib import Path
import numpy as np
import scipy
from scipy.optimize import minimize_scalar
from scipy.stats import qmc

A=3/8
G=math.sqrt(3)/2
W=-1/4+1j*math.sqrt(3)/4
D0=1/(4*math.sqrt(2))


def phase(z,d):
    y=z**4/4
    return 1j*z-y + d*(1j*z-2*y-W)**2/(2*(1+2*d*y))


def required(z,d):
    z4=z**4
    real=(-z4/4+d*(-z4/8+1/32-(z-G/2)**2/2))/(1+d*z4/2)
    return z*z*(A+real)


def exact_enough(action,d,grid):
    vals=action(grid,d)
    imax=int(np.argmax(vals))
    if 0<imax<len(grid)-1:
        ans=minimize_scalar(lambda x:-action(x,d),
                            bounds=(grid[imax-1],grid[imax+1]),method='bounded',
                            options={'xatol':2e-14})
        return float(ans.x),-float(ans.fun)
    return float(grid[imax]),float(vals[imax])


def audit(shard,repeat,power):
    u=qmc.Sobol(5,scramble=True,rng=np.random.default_rng(392689+shard+15*repeat)).random_base2(power)
    eps,din,p,s,q=u.T
    x=-3+6*p
    d=np.exp(np.log(2.**-16)+((shard+eps)/15)*np.log(.75/(2.**-16)))
    D=D0*(.82+.20*din)
    error=required(x,d)-D
    re_h=phase(x+.5j,d).real+A
    re_v1=phase(3+.5j*s,d).real+A
    re_v2=phase(-3+.5j*q,d).real+A
    direct=x*x*(A+phase(x,d).real)
    disagreement=float(np.max(np.abs(direct-required(x,d))))
    assert disagreement<1e-10
    return dict(points=len(u),error_positive=int(np.count_nonzero(error>0)),
                horizontal_positive=int(np.count_nonzero(re_h>1e-12)),
                max_horizontal=float(np.max(re_h)),max_vertical=float(np.max(np.maximum(re_v1,re_v2))),
                max_real_formula_difference=disagreement,
                sampled_range=[float(min(d)),float(max(d))])


def main():
    p=argparse.ArgumentParser()
    p.add_argument('shard',type=int)
    p.add_argument('--exponent',type=int,default=19)
    p.add_argument('--replicates',type=int,default=4)
    p.add_argument('--out',type=Path,required=True)
    args=p.parse_args()
    started=time.monotonic()
    cpu=time.process_time()
    rows=[audit(args.shard,i,args.exponent) for i in range(args.replicates)]
    dg=np.exp(np.log(2.**-16)+(args.shard+np.linspace(0,1,256))/15*np.log(.75/(2.**-16)))
    probe=np.linspace(0,3,1025)
    verts=np.linspace(-3,3,2049)
    probes=[]
    for d in dg:
        z,val=exact_enough(required,d,probe)
        nv=float(max(phase(verts+.5j,d).real)+A)
        probes.append(dict(delta=float(d),argmax=z,required=val,max_sampled_horizontal_excess=nv))
    best=min(probes,key=lambda r:r['required'])
    eligible=[r for r in probes if r['max_sampled_horizontal_excess']<1e-10]
    out=dict(shard=args.shard,pid=os.getpid(),points=sum(r['points'] for r in rows),rows=rows,
             best_error_threshold=best,best_horizontal_eligible=min(eligible,key=lambda r:r['required']) if eligible else None,
             tested_deltas=len(probes),entries=probes,elapsed=time.monotonic()-started,
             cpu=time.process_time()-cpu,python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)
    args.out.write_text(json.dumps(out,indent=2)+"\n")


if __name__=='__main__':
    main()
