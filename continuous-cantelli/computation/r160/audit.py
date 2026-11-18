#!/usr/bin/env python3
"""Finite-stage exact interval formula and comparisons with Metal output."""
import json
import math
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
import scipy
from scipy.optimize import brentq
from scipy.special import ndtr

NORM=1/np.sqrt(2*np.pi)
TARGET=np.sqrt(8/7)*np.exp(-2/7)


def geom(n):
    centers=np.array([0.])
    gaps=[]
    length=1.
    for k in range(1,n+1):
        w=length/(256*2**k)
        gaps.append((centers.copy(),w))
        new=(length-w)/2
        shift=(length-new)/2
        centers=np.c_[centers-shift,centers+shift].ravel()
        length=new
    return centers-length/2,centers+length/2,gaps,length


def responses(x,v,l,r):
    root=math.sqrt(v)
    a=(x-l)/root
    b=(x-r)/root
    p=NORM*np.exp(-a*a/2)/root
    q=NORM*np.exp(-b*b/2)/root
    return np.array([np.sum(ndtr(a)-ndtr(b)),sum(p-q),sum(b*q-a*p)/root])


def root_for(l,r,mult=1.):
    return brentq(lambda v:mult*responses(0,v,l,r)[0]-TARGET,.125,1/6,xtol=2e-15)


def central_gpu_root(queries,values):
    isort=np.argsort(queries[:,1])
    x=queries[isort,1]
    y=values[isort,0]
    hit=np.flatnonzero((y[:-1]>=TARGET)&(y[1:]<TARGET))
    assert len(hit)==1
    j=hit[0]
    return x[j]+(x[j+1]-x[j])*(y[j]-TARGET)/(y[j]-y[j+1])


def gap_value(gaps,x,v,order):
    total=0.
    for center,w in gaps:
        d=center-x
        if order=="midpoint":
            res=(NORM/np.sqrt(v))*w*np.exp(-(d*d)/(2*v))
        else:
            res=ndtr((d+w/2)/np.sqrt(v))-ndtr((d-w/2)/np.sqrt(v))
        total+=sum(res)
    return total


def main():
    tag=sys.argv[1] if len(sys.argv)>1 else "main"
    rec=json.loads(Path(f"receipt-{tag}.json").read_text())
    cert=json.loads(Path("certificate.json").read_text())
    vstar=float((Fraction(cert["v_lo"])+Fraction(cert["v_hi"]))/2)
    mass=np.prod(1.0-1.0/(256.0*np.exp2(np.arange(1,101))))
    Q=np.fromfile(f"queries-{tag}.f32",dtype="<f4").reshape((-1,4)).astype(float)
    V=np.fromfile(f"values-{tag}.f32",dtype="<f4").reshape((-1,6)).astype(float)
    assert len(Q)==len(V)==rec["rows"] and np.isfinite(V).all()
    assert np.isfinite(Q).all()
    coarse_n=rec["coarseX"]*rec["coarseV"]*len(rec["domainStages"])
    records=[]
    maximum_error=np.zeros(3)
    for n in rec["centerStages"]:
        l,r,gaps,width=geom(n)
        exact_root=root_for(l,r)
        adjusted_root=root_for(l,r,mass/(width*2**n))
        root_error=exact_root-vstar
        finer={}
        for order in (1,4):
            ids=np.flatnonzero((Q[coarse_n:,0]==0)&(Q[coarse_n:,2]==n)&(Q[coarse_n:,3]==order))+coarse_n
            center=central_gpu_root(Q[ids],V[ids])
            j=ids[np.argmin(abs(Q[ids,1]-vstar))]
            point=Q[j,1]
            real=responses(0,point,l,r)
            base=gap_value(gaps,0,point,"analytic")
            midpoint=gap_value(gaps,0,point,"midpoint")
            finer[f"order_{order}"]=dict(root=float(center),gpu_minus_exact_stage_root=float(center-exact_root),
                  checked_v=point,gpu_minus_cpu=tuple(V[j,:3]-real),gpu_gap_mass=V[j,3],
                  exact_gap_mass=float(base),midpoint_gap_error=float(midpoint-base),
                  gpu_minus_exact_gap=float(V[j,3]-base),latest_gap_mass=V[j,5])
            maximum_error=np.maximum(maximum_error,np.abs(V[j,:3]-real))
        space={}
        if n in rec["domainStages"]:
            ids=np.flatnonzero(Q[:coarse_n,2]==n)
            nv,nx=rec["coarseV"],rec["coarseX"]
            values=V[ids].reshape(nv,nx,6)
            query=Q[ids].reshape(nv,nx,4)
            argmax=query[np.arange(nv),np.argmax(values[:,:,0],axis=1),0]
            assert np.all(argmax==0)
            curvature=float(np.max(values[:,:,2][np.abs(query[:,:,0])<=.5]))
            assert curvature<0
            screen_v=[0,nv//2,nv-1]
            screen_x=[0,nx//4,nx//2,3*nx//4,nx-1]
            for i in screen_v:
                for j in screen_x:
                    x,v=query[i,j,:2]
                    real=responses(x,v,l,r)
                    maximum_error=np.maximum(maximum_error,np.abs(values[i,j,:3]-real))
            space=dict(variance_count=nv,sample_count=nx,all_sampled_argmax_zero=True,
                       maximum_curvature_on_core=curvature)
        records.append(dict(stage=n,root=exact_root,root_minus_infinite=root_error,
                            scaled_root_error=(2**n)*root_error,
                            mass_rescaled_root=adjusted_root,mass_rescaled_minus_infinite=adjusted_root-vstar,
                            cpu_stage_curvature=responses(0,exact_root,l,r)[2],center=finer,spatial=space))
    counts=dict(gap_kernel_nodes=int(np.sum(((1<<Q[:,2].astype(np.int64))-1)*Q[:,3])),
                base_kernel_nodes=16*len(Q))
    result=dict(ngrid=len(Q),kernel_node_counts=counts,max_abs_v_dx_dxx_error=maximum_error.tolist(),
                infinite_v=vstar,scipy=scipy.__version__,numpy=np.__version__,stages=records)
    Path(f"audit-{tag}.json").write_text(json.dumps(result,indent=2)+"\n")
    print("counts",counts)
    print("Float32 max absolutes [u,u_x,u_xx]",maximum_error)
    for r in records:
        print(r["stage"],r["root"],r["scaled_root_error"],r["mass_rescaled_minus_infinite"],
              r["center"]["order_4"]["gpu_minus_exact_stage_root"],
              r["center"]["order_4"]["midpoint_gap_error"],
              r["center"]["order_4"]["gpu_minus_exact_gap"],
              r["spatial"].get("maximum_curvature_on_core"))


if __name__=="__main__":
    main()
