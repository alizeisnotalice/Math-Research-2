import argparse
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as Q
import json
import multiprocessing
import os
from pathlib import Path
import platform
import re
import statistics
import subprocess
from time import monotonic, sleep

import bound_primitives as b

ROOT=Path(__file__).parent
S=b.SCALE
WW=b.WEIGHT_DEN
R=b.R
DEFAULT_A=('171/100','437/256','6999/4096')


def dist(ix,c):
    a,z=ix-c,ix+1-c
    return (min(abs(a),abs(z)) if a*z>0 else 0),max(abs(a),abs(z))


def positive_tail():
    upper=0
    for kind,t,y,w in b.TESTS:
        if w<=0:
            continue
        if kind=='t':
            upper+=2*w*b.efrac((R-abs(t))**2/2-t*t/2)[1]
        elif kind=='h':
            inv=b.sr_inv(t.numerator,t.denominator)[1]
            upper+=2*w*b.ceildiv(inv*b.efrac(Q(R*R,2))[1],S)
        else:
            upper+=b.ceildiv(2*w*b.H.numerator*b.efrac(Q(R*R,2))[1],b.H.denominator)
    return Q(upper,S*WW)


def worker(payload):
    shard,nx,nf,anchors=payload
    started=b.utc()
    mark=monotonic()
    lo=[0]*len(anchors)
    up=[0]*len(anchors)
    boxes=leaves=points=kept=0
    times=[(int(t),w,b.epositive(t*t/2)) for kind,t,_,w in b.TESTS if kind=='t']
    heats=[(int((b.H+s)*nf),int(y*nx),w) for kind,s,y,w in b.TESTS if kind=='h']
    roots=[None]+[b.sr_inv(i,nf) for i in range(1,max(v for v,_,_ in heats)+1)]
    cap=int(b.H*nf)
    xden=2*nx*nx
    fscale=nx*nx//nf
    exp_lo,exp_hi,ceildiv=b.exp_lower,b.exp_upper,b.ceildiv
    for ix in range(-R*nx+shard,R*nx,12):
        xmn,xmx=b.square_bounds(ix,ix+1,0)
        gl,gu=exp_lo(xmx,xden),exp_hi(xmn,xden)
        tv=[]
        for t,w,ep in times:
            mn,mx=b.square_bounds(ix,ix+1,t*nx)
            tv.append((t*t*fscale,w,ep,mn,mx))
        hv=[]
        for var,y,w in heats:
            mn,mx=b.square_bounds(ix,ix+1,y)
            hv.append((var,w,mn*nf,mx*nf))

        def lower(flo,fhi):
            nonlocal boxes
            boxes+=1
            total=WW*flo*gl//nf
            for tsq,w,ep,mn,mx in tv:
                if w>0:
                    val=exp_lo(mx+tsq*fhi,xden)*ep[0]//S
                else:
                    val=ceildiv(exp_hi(mn+tsq*flo,xden)*ep[1],S)
                total+=w*val
            for var,w,mn,mx in hv:
                vlo,vhi=var-fhi,var-flo
                if w>0:
                    val=exp_lo(mx,xden*vlo)*roots[vhi][0]*gl//(S*S)
                else:
                    val=ceildiv(exp_hi(mn,xden*vhi)*roots[vlo][1]*gu,S*S)
                total+=w*val
            return total

        def point(f):
            nonlocal points
            points+=1
            total=ceildiv(WW*f*gu,nf)
            for tsq,w,ep,mn,mx in tv:
                if w>0:
                    val=ceildiv(exp_hi(mn+tsq*f,xden)*ep[1],S)
                else:
                    val=exp_lo(mx+tsq*f,xden)*ep[0]//S
                total+=w*val
            for var,w,mn,mx in hv:
                v=var-f
                if w>0:
                    val=ceildiv(exp_hi(mn,xden*v)*roots[v][1]*gu,S*S)
                else:
                    val=exp_lo(mx,xden*v)*roots[v][0]*gl//(S*S)
                total+=w*val
            return total

        def bracket(start,feas):
            nonlocal leaves,kept
            # The incumbent comes from a point satisfying the envelope over
            # the whole spatial cell. A lower box above it can be pruned.
            samples=[feas+(cap-feas)*k//8 for k in range(9)]
            incumbent=min(point(p) for p in samples)
            stack=[(start,cap,lower(start,cap))]
            best_lower=None
            while stack:
                a,z,price=stack.pop()
                if price>=incumbent:
                    kept+=1
                    continue
                if z-a<=1:
                    leaves+=1
                    if z>=feas:
                        incumbent=min(incumbent,point(z))
                    if best_lower is None or price<best_lower:
                        best_lower=price
                    continue
                mid=(a+z)//2
                if mid>=feas:
                    incumbent=min(incumbent,point(mid))
                left=lower(a,mid)
                right=lower(mid,z)
                if left<=right:
                    stack.extend(((mid,z,right),(a,mid,left)))
                else:
                    stack.extend(((a,mid,left),(mid,z,right)))
            return min(best_lower,incumbent) if best_lower is not None else incumbent,incumbent

        local=(-3*nx//4<=ix<3*nx//4)
        if local:
            dmin,dmax=dist(ix,nx//2)
            for i,a in enumerate(anchors):
                num,den=a.numerator,a.denominator
                start=max(0,(num*nx-dmax*den)*nf//(den*nx))
                feas=max(0,ceildiv((num*nx-dmin*den)*nf,den*nx))
                v,w=bracket(start,feas)
                lo[i]+=v
                up[i]+=w
        else:
            v,w=bracket(0,0)
            for i in range(len(anchors)):
                lo[i]+=v
                up[i]+=w
    return dict(pid=os.getpid(),shard=shard,started=started,ended=b.utc(),seconds=monotonic()-mark,
                lower=lo,upper=up,boxes=boxes,finest_leaves=leaves,point_enclosures=points,pruned=kept)


def monitor(pids):
    ps=subprocess.run(['ps','-p',','.join(map(str,pids)),'-o','pid=,%cpu=,state=,etime='],
                      capture_output=True,text=True)
    stamp=b.utc()
    top=subprocess.run(['top','-l','2','-s','1','-n','0'],capture_output=True,text=True)
    return dict(utc=stamp,ps=ps.stdout.splitlines(),load=[s for s in top.stdout.splitlines()
                 if 'CPU usage:' in s][-1:])


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--nx',type=int,default=4096)
    parser.add_argument('--nf',type=int,default=4096)
    parser.add_argument('--prefix',default='target')
    parser.add_argument('--exp-grid-log2',type=int,default=12)
    parser.add_argument('anchor',nargs='*',default=DEFAULT_A)
    args=parser.parse_args()
    nx,nf=args.nx,args.nf
    assert nx%4==nf%4==0 and nx*nx%nf==0
    anchors=tuple(Q(a) for a in args.anchor)
    assert all(0<a<b.H for a in anchors)
    launch=b.utc()
    b.EXP_GRID=1<<args.exp_grid_log2
    b.EXP_END=40*b.EXP_GRID
    b.exp_table()
    with ProcessPoolExecutor(max_workers=12,mp_context=multiprocessing.get_context('fork')) as pool:
        jobs=[pool.submit(worker,(k,nx,nf,anchors)) for k in range(12)]
        processes=pool._processes
        pids=list(processes)
        meta=dict(start=launch,parent=os.getpid(),pids=pids,nx=nx,nf=nf,
                  exponential_grid=b.EXP_GRID,platform=platform.platform(),python=platform.python_version())
        tele=ROOT/f'{args.prefix}-telemetry.jsonl'
        with tele.open('w') as f:
            print(json.dumps(meta),file=f,flush=True)
            while not all(j.done() for j in jobs):
                print(json.dumps(monitor(pids)),file=f,flush=True)
                sleep(.25)
        workers=[j.result() for j in jobs]
    ps=subprocess.run(['ps','-p',','.join(map(str,pids)),'-o','pid='],capture_output=True,text=True)
    tlo,thi,tn=b.targets_and_tail()
    tp=positive_tail()
    gl,gu=b.gamma_normalizer()
    checks=[]
    for i,a in enumerate(anchors):
        rawlo=sum(w['lower'][i] for w in workers)
        rawup=sum(w['upper'][i] for w in workers)
        assert 0<rawlo<=rawup
        lo=Q(rawlo,nx*WW*gu)-Q(thi+tn,S*WW)
        hi=Q(rawup,nx*WW*gl)-Q(tlo,S*WW)+tp
        assert lo<=hi
        checks.append(dict(a=str(a),lower=str(lo),upper=str(hi),low=float(lo),high=float(hi)))
    samples=[json.loads(l) for l in tele.read_text().splitlines()[1:]]
    full=[m for m in samples if len(m['ps'])==12 and all(w.split()[2][0]=='R' for w in m['ps'])]
    group=[sum(float(w.split()[1]) for w in m['ps']) for m in full]
    idle=[float(re.search(r'([\d.]+)% idle',m['load'][0])[1]) for m in full]
    out=meta|dict(ended=b.utc(),workers=workers,brackets=checks,
                  normalization=[gl,gu],target_raw=[tlo,thi],
                  tail_negative=str(Q(tn,S*WW)),tail_positive=str(tp),
                  stats={k:sum(w[k] for w in workers) for k in
                          ('boxes','finest_leaves','point_enclosures','pruned')},
                  all_R_samples=len(full),group_cpu_median=statistics.median(group) if group else None,
                  host_idle=[min(idle),max(idle)] if idle else None,
                  exits={pid:p.exitcode for pid,p in processes.items()},processes_after=ps.stdout.splitlines())
    (ROOT/f'{args.prefix}.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(checks),flush=True)


if __name__=='__main__':
    main()
