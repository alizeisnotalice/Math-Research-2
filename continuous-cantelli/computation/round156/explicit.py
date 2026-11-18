from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from time import monotonic
import json
import os
import subprocess

import numpy as np

from run import initial, local_sink, mesh, phi, utc

ROOT=Path(__file__).parent


def run(cfg):
    name,eps,n,radius=cfg
    started=utc()
    tick=monotonic()
    lo=2.
    width=eps**1.5/100
    hi=lo+width
    x=mesh(lo,hi,radius,n)
    mid=(lo+hi)/2
    data=initial(x,lo,hi,eps)
    values=data[1:-1].copy()
    sink=local_sink(x,lo,hi,eps)
    left,right=np.diff(x)[:-1],np.diff(x)[1:]
    lower=1/(left*(left+right))
    upper=1/(right*(left+right))
    dt=.78/max(lower+upper)
    dl,du=dt*lower,dt*upper
    diag=1-dl-du
    forcing=dt*sink
    x=x[1:-1]
    offsets=[-2,0,.125,.25,.375,.5,.625,.75,.875,1,3]
    ids=np.array([np.argmin(abs(x-(lo+a*width))) for a in offsets])
    first=np.full(len(ids),-1,int)
    basis=data[np.argmin(abs(x-mid))+1]/phi(mid)
    bound=int(1.10*basis/dt)+1
    pre_contact=np.zeros_like(first)
    recovery=0
    for i in range(1,bound+1):
        nxt=diag*values-forcing
        nxt[1:]+=dl[1:]*values[:-1]
        nxt[:-1]+=du[:-1]*values[1:]
        np.maximum(nxt,0.,out=nxt)
        values=nxt
        probe=values[ids]
        current=probe<=0
        newly=(first<0)&current
        first[newly]=i
        recovery+=int(np.sum((first>=0)&(~current)))
        if i>=bound*.96 and np.all(first>=0):
            break
    out=dict(name=name,eps=eps,pid=os.getpid(),started=started,ended=utc(),seconds=monotonic()-tick,
             nodes=len(x)+2,dt=dt,cfl=float(max(dl+du)),steps=i,recovery=recovery,
             crossings=[dict(offset=a,x=float(x[j]),bracket=[(step-1)*dt,step*dt] if step>0 else None)
                        for a,step,j in zip(offsets,first,ids)])
    (ROOT/f'{name}.json').write_text(json.dumps(out,indent=2)+'\n')
    return out


def main():
    cases=[('explicit16',.125,16,.08),('explicit24',.125,24,.08),('explicit_eps16',.25,16,.14)]
    with ProcessPoolExecutor(max_workers=len(cases)) as pool:
        work=[pool.submit(run,cfg) for cfg in cases]
        process=pool._processes
        result=[job.result() for job in work]
    pids=list(process)
    ps=subprocess.run(['ps','-p',','.join(map(str,pids)),'-o','pid='],
                      capture_output=True,text=True)
    info=dict(parent=os.getpid(),pids=pids,exits={k:v.exitcode for k,v in process.items()},
              remaining=ps.stdout.splitlines(),runs=result)
    (ROOT/'summary-explicit.json').write_text(json.dumps(info,indent=2)+'\n')
    print(json.dumps(info),flush=True)


if __name__=='__main__':
    main()
