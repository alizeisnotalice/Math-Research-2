import argparse
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import json
from math import exp, sqrt
import os
from pathlib import Path
import platform
import subprocess
from time import monotonic, sleep
import warnings

import numpy as np
import scipy
from scipy.optimize import OptimizeWarning, linprog
from scipy.sparse import csr_matrix, hstack
from scipy.special import ndtr

ROOT=Path(__file__).parent
TESTS=(('moment',0.,0.),('time',-2.,0.),('time',-1.,0.),('time',1.,0.),
       ('time',2.,0.),('heat',.25,-1.5),('heat',.25,-.5),('heat',.25,.5),('heat',.25,1.5),
       ('heat',.75,-1.),('heat',.75,1.))
EDGES=np.linspace(-8,8,129)
X=(EDGES[:-1]+EDGES[1:])/2
MASS=np.diff(ndtr(EDGES))
NQ=17


def utc():
    return datetime.now(timezone.utc).isoformat()


def feature(x,f,h,tests):
    output=[]
    for kind,a,b in tests:
        if kind=='moment':
            value=f-1
        elif kind=='time':
            value=np.exp(a*x-.5*a*a*f)-1
        else:
            var=h+a-f
            target=exp(-b*b/(2*(h+a)))/sqrt(h+a)
            value=np.exp(-(x-b)**2/(2*var))/np.sqrt(var)-target
        output.append(value)
    return np.stack(output,axis=-1)


def enclosing(a0,a1,delta,mod,r,h):
    lo=np.maximum(0.,a0+delta-mod*np.abs(X-.5))
    hi=np.minimum(h,a1+mod*np.abs(X+.5))
    return np.where(np.abs(X)<=r,lo,0),np.where(np.abs(X)<=r,hi,h)


def solve(a0,a1,delta,mod,r,h,tests=TESTS):
    lo,hi=enclosing(a0,a1,delta,mod,r,h)
    if (lo>hi).any():
        return None
    f=lo[:,None]+(hi-lo)[:,None]*np.linspace(0.,1.,NQ)
    mat=feature(X[:,None],f,h,tests).reshape(-1,len(tests))
    link=csr_matrix((np.ones(len(X)*NQ),(np.arange(len(X)*NQ),np.repeat(np.arange(len(X)),NQ))),
                    shape=(len(X)*NQ,len(X)))
    constr=hstack([csr_matrix(-mat),link],format='csr')
    objective=np.concatenate([np.zeros(len(tests)),-MASS])
    with warnings.catch_warnings():
        warnings.simplefilter('ignore',OptimizeWarning)
        sol=linprog(objective,A_ub=constr,b_ub=np.zeros(len(X)*NQ),
                    bounds=[(-1,1)]*len(tests)+[(None,None)]*len(X),
                    method='highs',options={'threads':1})
    if not sol.success:
        raise RuntimeError(sol.message)
    return -sol.fun,sol.x[:len(tests)].tolist()


def worker(arg):
    k,limit=arg
    start=utc()
    initial=monotonic()
    rows=[]
    for j in range(420):
        i=12*j+k
        h=[1.25,1.5,1.75,2.,2.5][(i//12)%5]
        rad=[.75,1.,1.25,1.5,2.][(i//60)%5]
        mod=[.5,.75,1.,1.5,2.,2.5][(i*11+i//7)%6]
        a0=((i*19+i//5)%13)/16
        a1=a0+1/16
        delta=[.25,.375,.5,.625,.75,1.][(i*13+i//11)%6]
        if delta>mod or a0+delta>h:
            continue
        result=solve(a0,a1,delta,mod,rad,h)
        if result is not None:
            margin,coeff=result
            row=dict(i=i,a0=a0,a1=a1,delta=delta,lip=mod,radius=rad,H=h,margin=margin,coeff=coeff)
            if margin>.005:
                time=solve(a0,a1,delta,mod,rad,h,tuple(t for t in TESTS if t[0]!='heat'))[0]
                heat=solve(a0,a1,delta,mod,rad,h,tuple(t for t in TESTS if t[0]!='time'))[0]
                row|=dict(time_margin=time,heat_margin=heat,joint_gain=margin-max(time,heat))
                if min(time,heat)<.002 and max(time,heat)<.005:
                    row['strictly_joint']=True
            rows.append(row)
        if monotonic()-initial>=limit:
            break
    end=utc()
    (ROOT/f'shard-{k:02d}.json').write_text(json.dumps(rows,separators=(',',':'))+'\n')
    return dict(pid=os.getpid(),shard=k,started=start,ended=end,n=len(rows))


def monitor(pids):
    ps=subprocess.run(['ps','-p',','.join(map(str,pids)),'-o','pid=,%cpu=,state=,etime='],
                      capture_output=True,text=True)
    stamp=utc()
    top=subprocess.run(['top','-l','2','-s','1','-n','0'],capture_output=True,text=True)
    return dict(utc=stamp,workers=ps.stdout.splitlines(),load=[l for l in top.stdout.splitlines() if 'CPU usage:' in l][-1:])


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--duration',type=float,default=33)
    config=ap.parse_args()
    start=utc()
    with ProcessPoolExecutor(max_workers=12) as pool:
        jobs=[pool.submit(worker,(i,config.duration)) for i in range(12)]
        process=pool._processes
        pids=list(process)
        meta=dict(start=start,parent=os.getpid(),pids=pids,python=platform.python_version(),
                  scipy=scipy.__version__,numpy=np.__version__,hardware=platform.platform(),tests=TESTS)
        with (ROOT/'telemetry.jsonl').open('w') as f:
            print(json.dumps(meta),file=f,flush=True)
            while not all(job.done() for job in jobs):
                print(json.dumps(monitor(pids)),file=f,flush=True)
                sleep(.3)
        stats=[job.result() for job in jobs]
    rows=[row for path in ROOT.glob('shard-*.json') for row in json.loads(path.read_text())]
    ps=subprocess.run(['ps','-p',','.join(map(str,pids)),'-o','pid='],capture_output=True,text=True)
    out=meta|dict(end=utc(),stats=stats,n=len(rows),positive=sum(r['margin']>.005 for r in rows),
                  best_margin=sorted(rows,key=lambda r:-r['margin'])[:8],
                  best_joint_gain=sorted([r for r in rows if 'joint_gain' in r],key=lambda r:-r['joint_gain'])[:8],
                  strictly_joint=sorted([r for r in rows if r.get('strictly_joint')],key=lambda r:-r['margin'])[:12],
                  exits={p:proc.exitcode for p,proc in process.items()},after=ps.stdout.splitlines())
    (ROOT/'search-summary.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ['start','end','n','positive','after']}|{'strict_joint_count':len(out['strictly_joint'])}))


if __name__=='__main__':
    main()
