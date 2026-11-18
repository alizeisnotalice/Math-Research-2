import argparse
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from math import ceil, pi, sqrt
import json
import os
from pathlib import Path
import platform
import subprocess
from time import monotonic

import numpy as np
import scipy
from scipy.linalg import solve_banded
from scipy.special import ndtr, roots_legendre

ROOT=Path(__file__).parent


def utc():
    return datetime.now(timezone.utc).isoformat()


def phi(x):
    return np.exp(-.5*x*x)/sqrt(2*pi)


def mesh(lo,hi,radius,n,cap=.00022,growth=1.065):
    width=hi-lo
    step=width/n
    fine=np.linspace(lo-3*width,hi+3*width,7*n+1)
    outward=[0.]
    available=radius-3.5*width
    while outward[-1]<available:
        outward.append(min(available,outward[-1]+step))
        step=min(cap,growth*step)
    out=np.asarray(outward[1:])
    return np.concatenate([fine[0]-out[::-1],fine,fine[-1]+out])


def initial(x,lo,hi,eps):
    # K_e(r)=2*sqrt(e)*int_0^infty v phi(v+abs(r)/sqrt(e)) dv.
    zv,wv=roots_legendre(64)
    vv=7*(zv+1)
    wv=7*wv
    z,w=roots_legendre(14)
    packed=[]
    for locations in np.array_split(x,max(1,ceil(len(x)/64))):
        split=np.clip(locations,lo,hi)
        result=np.zeros(locations.shape)
        for a,b in [(np.full_like(split,lo),split),(split,np.full_like(split,hi))]:
            span=(b-a)/2
            y=(a+b)[:,None]/2+span[:,None]*z
            shifted=np.abs(locations[:,None]-y)/sqrt(eps)
            kern=2*sqrt(eps)*(phi(shifted[:,:,None]+vv[None,None,:])@(wv*vv))
            result+=span*(kern*phi(y)@w)
        packed.extend((1-eps/2)*result)
    return np.array(packed)


def local_sink(x,lo,hi,eps):
    # Split each dual cell at both layer endpoints and integrate mass.
    l=(x[:-2]+x[1:-1])/2
    r=(x[1:-1]+x[2:])/2
    a=np.clip(lo,l,r)
    b=np.clip(hi,l,r)
    outside=(ndtr(a)-ndtr(l))+(ndtr(r)-ndtr(b))
    inside=ndtr(b)-ndtr(a)
    return (outside+(eps/2)*inside)/(r-l)


def time_step(v,active,lower,diag,upper,rhs):
    band=np.empty((3,len(rhs)))
    band[0,0]=band[2,-1]=0.
    for it in range(50):
        band[0,1:]=np.where(active[:-1],0.,upper[:-1])
        band[1]=np.where(active,1.,diag)
        band[2,:-1]=np.where(active[1:],0.,lower[1:])
        v=solve_banded((1,1),band,np.where(active,0.,rhs),overwrite_ab=False,
                       overwrite_b=True,check_finite=False)
        residual=diag*v-rhs
        residual[1:]+=lower[1:]*v[:-1]
        residual[:-1]+=upper[:-1]*v[1:]
        new=(v<=residual)
        if np.array_equal(active,new):
            break
        active=new
    else:
        raise RuntimeError('nonconvergent contact policy')
    projected=np.maximum(v,0.)
    defect=max(float(np.abs(np.minimum(v,residual)).max()),
               max(0.,-float(v.min()),-float(residual.min())))
    return projected,active,it+1,defect


def solve(settings):
    label,eps,n,dt_ratio,radius,bcondition,*mesh_config=settings
    stamp=utc()
    tick=monotonic()
    lo=2.
    width=eps**1.5/100
    hi=lo+width
    x=mesh(lo,hi,radius,n,*mesh_config)
    d0=initial(x,lo,hi,eps)
    values=d0[1:-1].copy()
    sink=local_sink(x,lo,hi,eps)
    left=np.diff(x)[:-1]
    right=np.diff(x)[1:]
    ref=d0[np.argmin(abs(x-(lo+hi)/2))]/phi((lo+hi)/2)
    dt=ref*dt_ratio
    ldt=dt/((left+right)*left)
    rdt=dt/((left+right)*right)
    lower,diag,upper=-ldt,1+ldt+rdt,-rdt
    boundary=np.zeros_like(values)
    if bcondition=='frozen':
        boundary[0]=ldt[0]*d0[0]
        boundary[-1]=rdt[-1]*d0[-1]
    test_offsets=[-20,-10,-6,-3,-2,-1,-.5,-.25,0,.0625,.125,.25,.375,.5,
                  .625,.75,.875,.9375,1,1.25,1.5,2,3,4,7,11,21]
    actual=[int(np.argmin(np.abs(x[1:-1]-(lo+a*width)))) for a in test_offsets]
    near=np.abs(x[1:-1]-(lo+hi)/2)<=4*width
    target=set(actual)
    crossings=np.full(len(values),-1,dtype=int)
    active=np.zeros(len(values),dtype=bool)
    recoveries=0
    increases=0.
    max_iter=0
    count_iter=0
    max_defect=0.
    frames=[]
    for k in range(1,ceil(2.6/dt_ratio)+1):
        prev=values
        rhs=values-dt*sink+boundary
        values,active,it,defect=time_step(values,active,lower,diag,upper,rhs)
        max_iter=max(max_iter,it)
        count_iter+=it
        max_defect=max(max_defect,defect)
        increases=max(increases,float((values-prev).max()))
        recoveries+=int(np.sum((crossings>=0)&(~active)))
        new=(crossings<0)&active
        crossings[new]=k
        if k%max(1,int(.01/dt_ratio))==0:
            frames.append([float(k*dt),float(values[actual[len(actual)//2]]),
                          float(values[near].max()),int(active.sum())])
        if k*dt>=1.08*ref and all(crossings[j]>=0 for j in target):
            break
    position=x[1:-1]
    checks=[]
    for a,j in zip(test_offsets,actual):
        cross=int(crossings[j])
        bracket=[(cross-1)*dt,cross*dt] if cross>0 else None
        checks.append(dict(offset=a,x=float(position[j]),crossing_step=cross,bracket=bracket))
    near_pos=[[float(position[j]),int(crossings[j])] for j in np.flatnonzero(near)]
    finite=[j for j in np.flatnonzero(near) if crossings[j]>0]
    last=max(finite,key=lambda j:crossings[j]) if finite else None
    out=dict(label=label,eps=eps,points=len(x),K=[lo,hi],layer_width=width,
             inside_subdivisions=n,min_step=float(np.diff(x).min()),dt=dt,steps=k,reference_time=ref,
             radius=radius,boundary=bcondition,started=stamp,ended=utc(),seconds=monotonic()-tick,
             pid=os.getpid(),policy_iterations_max=max_iter,policy_iterations_mean=count_iter/k,
             lcp_defect_max=max_defect,max_positive_increase=increases,recovered_after_contact=recoveries,
             initial_at_test=[float(d0[j+1]) for j in actual],
             sink_at_test=[float(sink[j]) for j in actual],test_crossings=checks,frames=frames,
             termination_time=k*dt,
             near_last_contact=dict(x=float(position[last]),time=float(crossings[last]*dt)) if last is not None else None,
             near_crossings=near_pos)
    (ROOT/f'{label}.json').write_text(json.dumps(out,indent=2)+'\n')
    return {k:v for k,v in out.items() if k not in ('test_crossings','frames','near_crossings','initial_at_test','sink_at_test')}


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--mode',choices=['smoke','main','fine'],default='main')
    conf=p.parse_args()
    if conf.mode=='smoke':
        tasks=[('smoke',.125,16,.001,.06,'zero')]
    elif conf.mode=='fine':
        tasks=[('tiny_dt96',.125,96,1/86400,.08,'zero'),
               ('fine128',.125,128,1/57600,.08,'zero',.00011,1.045),
               ('fine192',.125,192,1/86400,.08,'zero',.00011,1.045),
               ('fine256',.125,256,1/172800,.1,'zero',.00011,1.04),
               ('wide_eps128',.25,128,1/48000,.14,'zero',.00011,1.045),
               ('wide_eps192',.25,192,1/96000,.14,'zero',.00011,1.045)]
    else:
        tasks=[('main16',.125,16,1/1800,.08,'zero'),('main32',.125,32,1/3600,.08,'zero'),
               ('main64',.125,64,1/7200,.08,'zero'),('main96',.125,96,1/10800,.08,'zero'),
               ('wide64',.125,64,1/7200,.12,'zero'),('narrow64',.125,64,1/7200,.045,'zero'),
               ('frozen64',.125,64,1/7200,.08,'frozen'),('refine_dt',.125,64,1/14400,.08,'zero'),
               ('wide_eps16',.25,16,1/2400,.14,'zero'),('wide_eps32',.25,32,1/4800,.14,'zero'),
               ('wide_eps64',.25,64,1/9600,.14,'zero'),('wide_eps_frozen',.25,64,1/9600,.14,'frozen')]
    with ProcessPoolExecutor(max_workers=len(tasks)) as pool:
        jobs=[pool.submit(solve,v) for v in tasks]
        processes=pool._processes
        meta=dict(parent=os.getpid(),pids=list(processes),start=utc(),python=platform.python_version(),
                  numpy=np.__version__,scipy=scipy.__version__,platform=platform.platform())
        print(json.dumps(meta),flush=True)
        rows=[j.result() for j in jobs]
    query=subprocess.run(['ps','-p',','.join(map(str,meta['pids'])),'-o','pid='],
                         capture_output=True,text=True)
    out=meta|dict(end=utc(),results=rows,exits={p:x.exitcode for p,x in processes.items()},
                  processes_after=query.stdout.splitlines())
    (ROOT/f'summary-{conf.mode}.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out),flush=True)


if __name__=='__main__':
    main()
