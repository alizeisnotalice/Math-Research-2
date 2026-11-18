import argparse
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from functools import lru_cache
import json
from math import ceil, exp, log, pi, sqrt
import os
from pathlib import Path
import platform
import subprocess
from time import monotonic, sleep

import numpy as np
import scipy
from scipy.integrate import quad
from scipy.special import gamma, iv, roots_legendre

ROOT=Path(__file__).parent
LOG2=log(2)
SQRT2=sqrt(2)
BETA=[.55,.8,1.,1.1,1.25,1.4,1.5,1.6,1.8,2.]


def utc():
    return datetime.now(timezone.utc).isoformat()


@lru_cache(maxsize=64)
def rule(lo,hi,npanel,nleg):
    z,w=roots_legendre(nleg)
    width=(hi-lo)/npanel
    nodes=lo+width*(np.arange(npanel)[:,None]+(z[None,:]+1)/2)
    return nodes.ravel(),np.tile((width/2)*w,npanel)


def geometry(log_t,beta,omega,theta,shift,amp,c):
    ell0=log_t-log(2*log_t)
    ell=ell0+shift
    relative=np.log1p(shift/ell0)
    phase=theta+omega*ell0**beta*np.expm1(beta*relative)
    bump=amp*np.cos(phase)
    offset=shift-log(2*log_t)+.5*np.log(ell)+c+bump
    return ell,log_t+offset,offset,bump


def point_integrand(w,log_t,beta,omega,theta,amp,c):
    y=2*log_t+sqrt(2*log_t)*w
    shift=-np.log1p(w/sqrt(2*log_t))
    _,d,offset,_=geometry(log_t,beta,omega,theta,shift,amp,c)
    return sqrt(log_t/pi)*np.exp(offset-(y-2*d)**2/(4*d) - .5*y*y*exp(-2*log_t))


def average_integrand(s,log_t,beta,omega,theta,amp,c):
    y=2*log_t*np.exp(s)
    ell,d,_,bump=geometry(log_t,beta,omega,theta,-s,amp,c)
    den=1+y*y/(2*d)
    return np.sqrt(2*ell/(pi*den))*np.exp(c+bump-d/den-.5*y*y*exp(-2*log_t))


def integrate_grid(fun,lo,hi,freq,profile,target_phase):
    n=ceil((hi-lo)*(1.4+0.80*freq))
    nodes,weights=rule(lo,hi,n,18)
    answer=float(weights@fun(nodes,*profile))
    if target_phase:
        nodes,weights=rule(lo,hi,2*n,22)
        finer=float(weights@fun(nodes,*profile))
        assert abs(answer-finer)<2e-7,(answer,finer,profile,n)
    return answer, n


def moments(log_t,beta,omega,theta,amp,c):
    ell0=log_t-log(2*log_t)
    outer_freq=omega*beta*ell0**(beta-1)
    point_freq=outer_freq/sqrt(2*log_t)
    n=np.arange(1,21)
    coeff=iv(n,amp)
    mix=iv(0,amp)
    exp_phase=np.exp(1j*n*theta)
    point_pred=exp(c)/SQRT2*(mix+2*np.sum(coeff*np.real(exp_phase)*np.exp(-.5*n*n*point_freq*point_freq)))
    # Linearization of the phase over the log-half-normal shift. Its validity
    # for beta>1 is not needed for the analytic convergence to the zero mode.
    freq=n*outer_freq
    mult=np.zeros(len(n),dtype=complex)
    good=freq<450
    mult[good]=np.exp(.5j*LOG2*freq[good])*gamma((1+1j*freq[good])/2)/sqrt(pi)
    avg_pred=exp(c)/SQRT2*(mix+2*np.sum(coeff*np.real(exp_phase*mult)))
    return point_freq,outer_freq,float(point_pred),float(avg_pred)


def evaluate(i,forced=None,check=None):
    # Each geometry is identified by the same center phase for both integrals.
    if forced is None:
        beta=BETA[(i//12)%len(BETA)]
        log_t=2**(6+((i*7+1)%33)/3)
        omega=.18+1.6*((i*sqrt(3)+.42)%1)
        theta=2*pi*((i*sqrt(5)+.14)%1)
        amp=.10+1.5*((i*sqrt(7)+.23)%1)
        ell0=log_t-log(2*log_t)
        omega=min(omega,55*sqrt(2*log_t)/(beta*ell0**(beta-1)))
    else:
        log_t,beta,omega,theta,amp=forced
    c=LOG2/2-log(iv(0,amp))  # Sets the common zero Fourier mode to one.
    point_freq,outer_freq,pred_p,pred_a=moments(log_t,beta,omega,theta,amp,c)
    params=(log_t,beta,omega,theta,amp,c)
    point,panels=integrate_grid(point_integrand,-9,11,point_freq,params,i%13==0)
    avg=None
    if outer_freq<=105:
        avg,n_avg=integrate_grid(average_integrand,-4,23,outer_freq,params,i%13==0)
    independent=None
    if check is True or (forced is None and i%179==0):
        val_p=quad(point_integrand,-9,11,args=params,epsabs=8e-11,epsrel=8e-11,
                   limit=max(96,ceil(12*point_freq)))[0]
        val_a=(quad(average_integrand,-4,23,args=params,epsabs=8e-11,epsrel=8e-11,
                    limit=max(96,ceil(20*outer_freq)))[0] if avg is not None else None)
        independent=[abs(point-val_p),abs(avg-val_a) if avg is not None else None]
        assert independent[0]<8e-8 and (val_a is None or independent[1]<8e-8), (i,independent)
    return dict(i=i,L=log_t,beta=beta,omega=omega,theta=theta,amp=amp,point=point,average=avg,
                predicted_point=pred_p,predicted_average=pred_a,point_freq=point_freq,
                outer_freq=outer_freq,independent=independent)


def worker(payload):
    k,duration=payload
    started=utc()
    start=monotonic()
    rows=[]
    for j in range(60000):
        row=evaluate(j*12+k)
        rows.append(row)
        if monotonic()-start>=duration:
            break
    path=ROOT/f'shard-{k:02d}.json'
    end=utc()
    path.write_text(json.dumps(rows,separators=(',',':'))+'\n')
    return dict(k=k,pid=os.getpid(),start=started,end=end,rows=len(rows))


def monitor(pids):
    ps=subprocess.run(['ps','-p',','.join(map(str,pids)),'-o','pid=,%cpu=,state=,etime='],
                      capture_output=True,text=True)
    stamp=utc()
    top=subprocess.run(['top','-l','2','-s','1','-n','0'],capture_output=True,text=True)
    lines=[line for line in top.stdout.splitlines() if 'CPU usage:' in line]
    return dict(utc=stamp,ps=ps.stdout.splitlines(),load=lines[-1:])


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--duration',type=float,default=36)
    args=parser.parse_args()
    launch=utc()
    with ProcessPoolExecutor(max_workers=12) as pool:
        futures=[pool.submit(worker,(k,args.duration)) for k in range(12)]
        processes=pool._processes
        pids=list(processes)
        capture=dict(start=launch,parent=os.getpid(),pids=pids,python=platform.python_version(),
                     platform=platform.platform(),numpy=np.__version__,scipy=scipy.__version__)
        with (ROOT/'telemetry.jsonl').open('w') as stream:
            print(json.dumps(capture),file=stream,flush=True)
            while not all(f.done() for f in futures):
                print(json.dumps(monitor(pids)),file=stream,flush=True)
                sleep(.3)
        job=[f.result() for f in futures]
    after=subprocess.run(['ps','-p',','.join(map(str,pids)),'-o','pid='],
                         capture_output=True,text=True)
    records=[r for path in ROOT.glob('shard-*.json') for r in json.loads(path.read_text())]
    checks=[v for row in records if row['independent'] for v in row['independent'] if v is not None]
    report=capture|dict(finished=utc(),workers=job,n=len(records),n_avg=sum(r['average'] is not None for r in records),
                        n_independent=len(checks),maximum_independent_disagreement=max(checks),
                        exits={pid:p.exitcode for pid,p in processes.items()},workers_after=after.stdout.splitlines())
    (ROOT/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)


if __name__=='__main__':
    main()
