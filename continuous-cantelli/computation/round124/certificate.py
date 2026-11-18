from array import array
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from fractions import Fraction as Q
from math import factorial, isqrt
import argparse
import json
import multiprocessing
import os
from pathlib import Path
import subprocess
from time import monotonic, sleep

ROOT=Path(__file__).parent
SCALE=1<<60
EXP_GRID=1<<12
EXP_END=40*EXP_GRID
WEIGHT_DEN=1000000
H=Q(5,2)
A0,A1=Q(3,4),Q(13,16)
DELTA,MOD,RAD=Q(1),Q(1),Q(3,4)
R=8
TESTS=(('m',Q(0),Q(0),1000000),
       ('t',Q(-2),Q(0),34939),('t',Q(-1),Q(0),361498),
       ('t',Q(1),Q(0),-15824),('t',Q(2),Q(0),128578),
       ('h',Q(1,4),Q(-3,2),-588210),('h',Q(1,4),Q(-1,2),-1000000),
       ('h',Q(1,4),Q(1,2),1000000),('h',Q(1,4),Q(3,2),-964943),
       ('h',Q(3,4),Q(-1),-182022),('h',Q(3,4),Q(1),1000000))
X1,X2=Q(-1,2),Q(1,2)
LOWER=UPPER=None
RIGHT_ONLY=False


def utc():
    return datetime.now(timezone.utc).isoformat()


def ceildiv(a,b):
    return -(-a//b)


def exp_table():
    """Directed dyadic intervals for exp(-k/EXP_GRID)."""
    global LOWER,UPPER
    u=Q(1,EXP_GRID)
    odd=sum((-u)**j/factorial(j) for j in range(8))
    even=odd+u**8/factorial(8)
    lo=(odd.numerator*SCALE)//odd.denominator
    hi=ceildiv(even.numerator*SCALE,even.denominator)
    LOWER,UPPER=array('Q',[SCALE]),array('Q',[SCALE])
    for _ in range(EXP_END):
        LOWER.append(LOWER[-1]*lo//SCALE)
        UPPER.append(ceildiv(UPPER[-1]*hi,SCALE))


def exp_lower(qnum,qden):
    i=ceildiv(qnum*EXP_GRID,qden)
    return LOWER[i] if i<=EXP_END else 0


def exp_upper(qnum,qden):
    i=qnum*EXP_GRID//qden
    return UPPER[min(i,EXP_END)]


def efrac(q):
    return exp_lower(q.numerator,q.denominator),exp_upper(q.numerator,q.denominator)


def epositive(q):
    lo,hi=efrac(q)
    return SCALE*SCALE//hi,ceildiv(SCALE*SCALE,lo)


def sr_inv(num,den):
    lower=isqrt(SCALE*SCALE*den//num)
    return lower,lower+1


def square_bounds(x0,x1,center):
    a,b=x0-center,x1-center
    return (min(a*a,b*b) if a*b>0 else 0),max(a*a,b*b)


def pi_bounds():
    def atan(z):
        odd=sum((-1)**j*z**(2*j+1)/(2*j+1) for j in range(22))
        return odd,odd+z**45/45
    xlo,xhi=atan(Q(1,5))
    ylo,yhi=atan(Q(1,239))
    return 16*xlo-4*yhi,16*xhi-4*ylo


def gamma_normalizer():
    lo,hi=pi_bounds()
    lo*=2
    hi*=2
    lower=isqrt(SCALE*SCALE*lo.numerator//lo.denominator)
    upper=isqrt(SCALE*SCALE*hi.numerator//hi.denominator)+1
    return lower,upper


def variance_indices(nx,nf,ix):
    x0,x1=Q(ix,nx),Q(ix+1,nx)
    if x0 < -RAD or x1 > RAD:
        return 0,int(H*nf)
    d1=max(abs(x0-X1),abs(x1-X1))
    d2=max(abs(x0-X2),abs(x1-X2))
    if RIGHT_ONLY:
        fmin=max(0,A0+DELTA-MOD*d2)
        fmax=H
    else:
        fmin=max(0,A0-MOD*d1,A0+DELTA-MOD*d2)
        fmax=min(H,A1+MOD*d1,A1+MOD*(X2-X1)+MOD*d2)
    fmin*=nf
    fmax*=nf
    return int(fmin),ceildiv(fmax.numerator,fmax.denominator)


def worker(payload):
    shard,nx,nf=payload
    start=utc()
    times=[(int(t),w,epositive(t*t/2)) for kind,t,_,w in TESTS if kind=='t']
    heat=[(int((H+s)*nf),int(y*nx),w) for kind,s,y,w in TESTS if kind=='h']
    sqr=[None]+[sr_inv(i,nf) for i in range(1,max(v for v,_,_ in heat)+1)]
    axden=2*nx*nx
    assert nx*nx%nf==0
    time_f_scale=nx*nx//nf
    h2den=2*nx*nx
    num=0
    boxes=0
    for ix in range(-R*nx+shard,R*nx,12):
        mn,mx=square_bounds(ix,ix+1,0)
        gl=exp_lower(mx,axden)
        gu=exp_upper(mn,axden)
        vtime=[(t*t*time_f_scale,w,ep,
                square_bounds(ix,ix+1,t*nx)[1 if w>0 else 0])
               for t,w,ep in times]
        vheat=[(var,w,square_bounds(ix,ix+1,y)[1 if w>0 else 0]*nf)
               for var,y,w in heat]
        f_lo,f_hi=variance_indices(nx,nf,ix)
        assert 0<=f_lo<f_hi<=int(H*nf)
        cell_min=None
        for flo in range(f_lo,f_hi):
            fhi=flo+1
            total=WEIGHT_DEN*flo*gl//nf
            for tsq,w,ep,dx2 in vtime:
                if w>0:
                    val=exp_lower(dx2+tsq*fhi,axden)*ep[0]//SCALE
                else:
                    val=ceildiv(exp_upper(dx2+tsq*flo,axden)*ep[1],SCALE)
                total+=w*val
            for var,w,dist2 in vheat:
                vlo,vhi=var-fhi,var-flo
                if w>0:
                    val=exp_lower(dist2,h2den*vlo)*sqr[vhi][0]*gl//(SCALE*SCALE)
                else:
                    val=ceildiv(exp_upper(dist2,h2den*vhi)*sqr[vlo][1]*gu,SCALE*SCALE)
                total+=w*val
            if cell_min is None or total<cell_min:
                cell_min=total
            boxes+=1
        num+=cell_min
    return dict(shard=shard,boxes=boxes,integral_numerator=num,pid=os.getpid(),
                start_utc=start,end_utc=utc())


def targets_and_tail():
    target_lo=target_hi=tail=0
    for kind,s,y,w in TESTS:
        if kind in ('m','t'):
            tlo=thi=SCALE
        else:
            v=H+s
            root=sr_inv(v.numerator,v.denominator)
            elo,ehi=efrac(y*y/(2*v))
            tlo=root[0]*elo//SCALE
            thi=ceildiv(root[1]*ehi,SCALE)
        if w>=0:
            target_lo+=w*tlo
            target_hi+=w*thi
        else:
            target_lo+=w*thi
            target_hi+=w*tlo
        if w>=0:
            continue
        if kind=='t':
            _,upper=efrac((R-abs(s))**2/2-s*s/2)
            tail+=2*abs(w)*upper
        elif kind=='h':
            _,sq=sr_inv(s.numerator,s.denominator)
            _,upper=efrac(Q(R*R,2))
            tail+=2*abs(w)*ceildiv(sq*upper,SCALE)
        else:
            _,upper=efrac(Q(R*R,2))
            tail+=ceildiv(2*abs(w)*H.numerator*upper,H.denominator)
    return target_lo,target_hi,tail


def mon(pids):
    query=subprocess.run(['ps','-p',','.join(map(str,pids)),'-o','pid=,%cpu=,state=,etime='],
                         capture_output=True,text=True)
    top=subprocess.run(['top','-l','2','-s','1','-n','0'],capture_output=True,text=True)
    return dict(utc=utc(),workers=query.stdout.splitlines(),load=[line for line in top.stdout.splitlines()
                                 if 'CPU usage:' in line][-1:])


def main():
    global RIGHT_ONLY
    ap=argparse.ArgumentParser()
    ap.add_argument('--nx',type=int,default=256)
    ap.add_argument('--nf',type=int,default=128)
    ap.add_argument('--right-anchor-only',action='store_true')
    ap.add_argument('--output-prefix',default='certificate')
    cfg=ap.parse_args()
    RIGHT_ONLY=cfg.right_anchor_only
    assert cfg.nx%4==cfg.nf%4==0 and cfg.nx*cfg.nx%cfg.nf==0
    initial=utc()
    exp_table()
    start=monotonic()
    with ProcessPoolExecutor(max_workers=12,mp_context=multiprocessing.get_context('fork')) as pool:
        futures=[pool.submit(worker,(i,cfg.nx,cfg.nf)) for i in range(12)]
        procs=pool._processes
        pids=list(procs)
        with (ROOT/f'{cfg.output_prefix}-telemetry.jsonl').open('w') as f:
            print(json.dumps(dict(start=initial,parent=os.getpid(),pids=pids)),file=f,flush=True)
            while not all(job.done() for job in futures):
                print(json.dumps(mon(pids)),file=f,flush=True)
                sleep(.2)
        results=[future.result() for future in futures]
    ps=subprocess.run(['ps','-p',','.join(map(str,pids)),'-o','pid='],
                      capture_output=True,text=True)
    lower_num=sum(w['integral_numerator'] for w in results)
    normal_lo,normal_hi=gamma_normalizer()
    assert lower_num>0
    integral_lower=Q(lower_num,cfg.nx*WEIGHT_DEN*normal_hi)
    target_lo,target_hi,tail=targets_and_tail()
    target=Q(target_hi,SCALE*WEIGHT_DEN)
    tail_bound=Q(tail,SCALE*WEIGHT_DEN)
    margin=integral_lower-target-tail_bound
    max_residual=margin/Q(sum(abs(w) for *_,w in TESTS),WEIGHT_DEN)
    out=dict(hardware='Apple M6',start=initial,end=utc(),duration=monotonic()-start,
             X_radius=R,nx=cfg.nx,nf=cfg.nf,boxes=sum(r['boxes'] for r in results),
             H=str(H),radius=str(RAD),lipschitz=str(MOD),anchor_points=[str(X1),str(X2)],
             left_anchor_interval=None if RIGHT_ONLY else [str(A0),str(A1)],
             required_increase=None if RIGHT_ONLY else str(DELTA),right_anchor_lower=str(A0+DELTA),
             tests=[[kind,str(s),str(y),str(Q(w,WEIGHT_DEN))] for kind,s,y,w in TESTS],
             raw=dict(integral_num=lower_num,scale=SCALE,weight_den=WEIGHT_DEN,
                      normal_lo=normal_lo,normal_hi=normal_hi,target_lo=target_lo,
                      target_hi=target_hi,tail=tail),
             certified_margin=str(margin),max_residual_lower=str(max_residual),
             numerical_readout=dict(integral_lower=float(integral_lower),target_upper=float(target),
                                    omitted_tail_upper=float(tail_bound),margin=float(margin),
                                    max_residual_lower=float(max_residual)),
             workers=results,exits={pid:p.exitcode for pid,p in procs.items()},remaining=ps.stdout.splitlines())
    (ROOT/f'{cfg.output_prefix}.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out['numerical_readout']),flush=True)
    assert margin>0


if __name__=='__main__':
    main()
