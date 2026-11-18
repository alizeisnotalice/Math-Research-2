import argparse
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import json
from math import ceil, exp, log, pi, sqrt
import os
from pathlib import Path
import platform
import subprocess
from time import monotonic, sleep

import numpy as np
import scipy
from scipy.optimize import brentq
from scipy.special import ndtr

ROOT = Path(__file__).parent
PHI = 1/sqrt(2*pi)
ROOT_TOL = 2e-13


def utc():
    return datetime.now(timezone.utc).isoformat()


def phi(x):
    return PHI*np.exp(-.5*x*x)


class Image:
    def __init__(self,a,b,t):
        self.a,self.b,self.sigma = a,b,sqrt(t)
        width=a+b
        nmax=ceil(5*self.sigma/width)+3
        centers=-2*width*np.arange(-nmax,nmax+1,dtype=float)
        self.centers=np.concatenate([centers,centers-2*a])
        self.sign=np.concatenate([np.ones_like(centers),-np.ones_like(centers)])

    def moment(self,left,right,order):
        s,c=self.sigma,self.centers
        lo,hi=(left-c)/s,(right-c)/s
        mass=np.where(lo>0,ndtr(-lo)-ndtr(-hi),ndtr(hi)-ndtr(lo))
        if order==0:
            return float(self.sign@mass)
        plo,phi_=phi(lo),phi(hi)
        first=c*mass+s*(plo-phi_)
        if order==1:
            return float(self.sign@first)
        second=c*c*mass+2*c*s*(plo-phi_)+s*s*(mass+lo*plo-hi*phi_)
        return np.array([self.sign@mass,self.sign@first,self.sign@second])


class Eigen:
    def __init__(self,a,b,t):
        self.a,self.b=a,b
        width=a+b
        nmax=ceil(sqrt(100*width*width/(pi*pi*t)))+4
        self.wave=np.arange(1,nmax+1,dtype=float)*pi/width
        self.factor=2/width*np.sin(self.wave*a)*np.exp(-.5*t*self.wave**2)

    def primitive(self,x,order):
        k=self.wave
        sine,cosine=np.sin(k*(x+self.a)),np.cos(k*(x+self.a))
        mass=-cosine/k
        if order==0:
            return mass
        first=x*mass+sine/k**2
        if order==1:
            return first
        second=x*x*mass+2*x*sine/k**2+2*cosine/k**3
        return np.array([mass,first,second])

    def moment(self,left,right,order):
        return (self.primitive(right,order)-self.primitive(left,order))@self.factor


def cap_opt(kern,mass,mean):
    a,b=kern.a,kern.b
    capmass=kern.moment(-a,b,0)
    eps=64*np.finfo(float).eps*max(1.,mass)
    if capmass < mass-eps:
        return dict(reason='mass',capmass=capmass)
    if capmass <= mass+eps:
        return dict(reason='near-mass-edge',capmass=capmass)
    cdf=lambda x: kern.moment(-a,x,0)
    def quantile(q):
        if q<=eps:
            return -a
        if q>=capmass-eps:
            return b
        return brentq(lambda x:cdf(x)-q,-a,b,xtol=ROOT_TOL)
    leftcut,rightcut=quantile(mass),quantile(capmass-mass)
    low,high=kern.moment(-a,leftcut,1),kern.moment(rightcut,b,1)
    if not low - eps <= mean <= high + eps:
        return dict(reason='mean',capmass=capmass,low=low,high=high)
    def cuts(q):
        return quantile(q),quantile(capmass-mass+q)

    if abs(mean-high)<=eps:
        q=0.
    elif abs(mean-low)<=eps:
        q=mass
    else:
        def residual(q):
            l,u=cuts(q)
            return kern.moment(-a,l,1)+kern.moment(u,b,1)-mean

        q=brentq(residual,0.,mass,xtol=ROOT_TOL)
    l,u=cuts(q)
    cap=kern.moment(-a,l,2)+kern.moment(u,b,2)
    gain=a*b*mass+(b-a)*mean-cap[2]
    return dict(reason='feasible',capmass=capmass,low=low,high=high,
                gain=float(gain),cuts=[l,u],mismatch=[float(cap[0]-mass),float(cap[1]-mean)])


def parameters(i):
    # Three irrational increments provide a repeatable low-discrepancy covering.
    a=exp(log(.24)+(log(3.6)-log(.24))*((.513+sqrt(2)*i)%1))
    b=exp(log(.24)+(log(3.6)-log(.24))*((.121+sqrt(3)*i)%1))
    time_ratio=exp(log(.002)+(log(.36)-log(.002))*((.882+sqrt(5)*i)%1))
    return a,b,time_ratio*(a+b)**2


def compare(i):
    a,b,t=parameters(i)
    mass=float(ndtr(-a)+ndtr(-b))
    mean=float(phi(b)-phi(a))
    base=float(a*(phi(b)-b*ndtr(-b))+b*(phi(a)-a*ndtr(-a)))
    e=cap_opt(Eigen(a,b,t),mass,mean)
    m=cap_opt(Image(a,b,t),mass,mean)
    assert e['reason']==m['reason'], (i,e,m)
    assert abs(m['capmass']-e['capmass'])<2e-12,(i,m,e)
    if e['reason']=='feasible':
        assert abs(m['gain']-e['gain'])<8e-10,(i,m,e)
        assert max(map(abs,m['mismatch']+e['mismatch']))<3e-11,(i,m,e)
        assert e['gain']>-1e-13 and e['gain']<=a*b*mass+(b-a)*mean+1e-10,(i,e)
        e['relative_gain']=e['gain']/base
        e['kernel_disagreement']=abs(m['gain']-e['gain'])
    return dict(i=i,a=a,b=b,t=t,mass=mass,mean=mean,base=base,**e)


def worker(args):
    shard,duration=args
    start=utc()
    timer=monotonic()
    result=[]
    for counter in range(100000):
        result.append(compare(12*counter+shard))
        if monotonic()-timer>=duration:
            break
    out=dict(shard=shard,pid=os.getpid(),start_utc=start,end_utc=utc(),rows=result)
    (ROOT/f'shard-{shard:02d}.json').write_text(json.dumps(out)+'\n')
    return {k:v for k,v in out.items() if k!='rows'}|{'n':len(result)}


def monitor(pids):
    query=subprocess.run(['ps','-p',','.join(map(str,pids)),'-o','pid=,%cpu=,state=,etime='],
                         capture_output=True,text=True)
    ps_time=utc()
    top=subprocess.run(['top','-l','2','-s','1','-n','0'],capture_output=True,text=True)
    load=[l for l in top.stdout.splitlines() if 'CPU usage:' in l]
    return dict(utc=ps_time,workers=query.stdout.splitlines(),load=load[-1:])


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--duration',type=float,default=30)
    args=parser.parse_args()
    start=utc()
    with ProcessPoolExecutor(max_workers=12) as pool:
        futures=[pool.submit(worker,(k,args.duration)) for k in range(12)]
        processes=pool._processes
        pids=list(processes)
        with (ROOT/'telemetry.jsonl').open('w') as stream:
            capture=dict(start_utc=start,parent=os.getpid(),pids=pids,platform=platform.platform(),
                         python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)
            print(json.dumps(capture),file=stream,flush=True)
            while not all(f.done() for f in futures):
                print(json.dumps(monitor(pids)),file=stream,flush=True)
                sleep(1.)
        shards=[f.result() for f in futures]
    query=subprocess.run(['ps','-p',','.join(map(str,pids)),'-o','pid='],capture_output=True,text=True)
    rows=[row for path in ROOT.glob('shard-*.json') for row in json.loads(path.read_text())['rows']]
    feasible=[row for row in rows if row['reason']=='feasible']
    answer=dict(**capture,end_utc=utc(),shards=shards,n=len(rows),feasible=len(feasible),
                mass_infeasible=sum(r['reason']=='mass' for r in rows),
                mean_infeasible=sum(r['reason']=='mean' for r in rows),
                near_mass_edge=sum(r['reason']=='near-mass-edge' for r in rows),
                kernel_disagreement=max(r['kernel_disagreement'] for r in feasible),
                relative_gain_quantile=np.quantile([r['relative_gain'] for r in feasible],
                                                   [0,.1,.5,.9,1]).tolist(),
                best=max(feasible,key=lambda r:r['relative_gain']),
                worker_exit_codes={pid:p.exitcode for pid,p in processes.items()},
                workers_after=query.stdout.splitlines())
    if len(feasible):
        (ROOT/'summary.json').write_text(json.dumps(answer,indent=2)+'\n')
    print(json.dumps(answer),flush=True)


if __name__=='__main__':
    main()
