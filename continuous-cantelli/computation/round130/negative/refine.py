"""New grids around the first negative real-space Bernstein test."""
import argparse
from datetime import datetime, UTC
import json, os, platform, time
from pathlib import Path
import numpy as np
import scipy
from scipy.optimize import minimize_scalar
from scan import grid, stable_results


GATE = datetime(2026, 11, 18, 7, 12, tzinfo=UTC).timestamp()


def shard(i, duration, gate, output):
    rng = np.random.default_rng(1211201+i)
    now = time.time()
    if now < gate:
        time.sleep(gate-now)
    start_time = time.time()
    start_mono = time.perf_counter()
    cpu = time.process_time()
    profiles = 0
    batches = 0
    best_list = []
    first_uncertain = None
    focus = (2.78, 3.10)
    wide = (-7.5, 9.5)
    full_j = np.arange(470+i, 640, 15)
    other_j = np.arange(1+i, 470, 15)
    top_samples = {}
    while time.perf_counter()-start_mono < duration or not batches:
        sample = batches % 5
        if sample == 0:
            js, interval, offset, count, order = full_j, wide, .5, 4608, 80
        elif sample == 1:
            js, interval, offset, count, order = other_j, wide, .5, 2304, 80
        else:
            js, interval, offset, count, order = full_j, focus, .65, 4096, (96, 128, 160)[sample-2]
        dx = (interval[1]-interval[0])/count
        y = interval[0]+(np.arange(count)+rng.uniform(0.,1.,count))*dx
        times, w, values = grid(y, order, 'three')
        for j in js:
            lo = max(0, int(offset*j))
            hi = min(1060-j, int((1.14 if sample<2 else .84)*j)+1)
            if j < 64:
                lo=0
                hi=min(1060-j,70)
            ks = np.arange(lo,hi+1)
            if not len(ks):
                continue
            raw, _ = stable_results(j,ks,times,w,values)
            profiles += len(ks)*len(y)
            if not np.any(raw < -1e-13):
                continue
            cols=raw.argmin(axis=1)
            values_at_min = raw[np.arange(len(ks)),cols]
            selected = np.where(values_at_min < -1e-13)[0][:2]
            for ki in selected:
                k=int(ks[ki]); yc=float(y[cols[ki]])
                cand=(int(j)+k, int(j), k, yc, float(values_at_min[ki]))
                previous=top_samples.get((j,k))
                if previous is None or cand[-1]<previous[-1]:
                    top_samples[j,k]=cand
            if time.perf_counter()-start_mono >= duration:
                break
        batches+=1
    end_time = time.time()
    end_cpu = time.process_time()-cpu
    for j,k in sorted(top_samples,key=lambda p:(sum(p),p[0]))[:30]:
        candidate=top_samples[j,k]
        def obj(y):
            t,w,v=grid(np.array([y]),192,'three')
            z,_=stable_results(j,np.array([k]),t,w,v)
            return float(z[0,0])
        yc=candidate[3]
        sol=minimize_scalar(obj,bounds=(yc-.001,yc+.001),method='bounded',options={'xatol':1e-11})
        best_list.append([int(j)+int(k),int(j),int(k),float(sol.x),float(sol.fun),obj(yc),candidate])
    rec=dict(worker=i,pid=os.getpid(),start_epoch=start_time,end_epoch=end_time,wall=end_time-start_time,
             cpu=end_cpu,seed=1211201+i,profiles=profiles,batches=batches,
             candidates=best_list,software=[platform.python_version(),np.__version__,scipy.__version__])
    output.write_text(json.dumps(rec,indent=2)+'\n')


def main():
    p=argparse.ArgumentParser()
    p.add_argument('worker',type=int)
    p.add_argument('out',type=Path)
    p.add_argument('--gate',type=float,default=GATE)
    p.add_argument('--duration',type=float,default=91.)
    args=p.parse_args()
    shard(args.worker,args.duration,args.gate,args.out)

if __name__=='__main__':
    main()
