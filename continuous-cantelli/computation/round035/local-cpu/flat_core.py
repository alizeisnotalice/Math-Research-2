"""Density and Taylor jets at a plateau/positive infinitely flat seam."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import json, math, os, platform, re, subprocess, sys, time
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.special import ndtr
from scipy.stats import qmc

R=0.25
SEED=301108
JETS=7
DPS=49
DIVISIONS=(2,3,4,5,6,7)
CPU_RE=re.compile(r"CPU usage: (.*?)% user, (.*?)% sys, (.*?)% idle")


def utc():
    return datetime.now(timezone.utc).isoformat(timespec="milliseconds")


def grid(size):
    rows=qmc.Sobol(6,scramble=True,seed=SEED).random_base2((size-1).bit_length())[:size]
    for j,(a,p,q,l,qq,ll) in enumerate(rows):
        yield j, float(4*a-.5), float(.7+1.7*p), float(.18+.3*q), float(2.8+3.7*l), \
            float(.18+.3*qq) if j%2 else 0.0, float(2.8+3.7*ll)


def hermite(k,z):
    a,b=1,z
    for j in range(1,k):
        a,b=b,z*b-j*a
    return a if k==0 else b


def panels(a,b,extras):
    return [a]+sorted({x for x in extras if a<x<b})+[b]


def baseline(a,y,e,side=None):
    if not e:
        return mp.exp(-a*y-y*y/2) if side is None or side*y>0 else mp.mpf(0)
    mean=(y-a*e)/(1+e)
    sd=mp.sqrt(e/(1+e))
    lo,hi=(-R,R) if side is None else ((0,R) if side>0 else (-R,0))
    def Phi(z):
        return mp.erfc(-z/mp.sqrt(2))/2
    density=mp.exp(a*a/2-(a+y)**2/(2*(1+e)))/mp.sqrt(1+e)
    return density*(Phi((hi-mean)/sd)-Phi((lo-mean)/sd))


def positive_mp(a,y,e,p,q,l,side,k=0):
    q0=1/mp.sqrt(2*mp.pi)
    yloc=side*y
    B=(q*R)**2
    le=mp.log(e) if e else None
    lB=mp.log(B)

    def logvar(u):
        return lB-l*((R/u)**p-1) if u else mp.ninf

    def f(u):
        if u<=0:
            return mp.mpf(0)
        lv=logvar(u)
        if e:
            lv=mp.log(e+mp.exp(lv)) if lv>le-120 else le
        dist= u-yloc
        if not dist:
            arg=mp.mpf(0)
        else:
            logarg=2*mp.log(abs(dist))-lv
            if logarg>20:
                return mp.mpf(0)
            arg=mp.exp(logarg)/2
        logkernel=-(k+1)*lv/2-arg-side*a*u-u*u/2
        if logkernel < -180:
            return mp.mpf(0)
        z=mp.sign(dist)*mp.sqrt(2*arg)
        return q0*mp.exp(logkernel)*hermite(k,z)*side**k

    ext=(R/16,R/8,R/4,R/2,R*3/4,R*7/8)
    if not 0<yloc<R:
        extra = ((3*abs(y),abs(y)) if y else ()) + (tuple(j*mp.sqrt(e) for j in (1,4,12)) if e else ())
        return mp.quad(f,panels(0,R,ext+extra))
    # Separate a narrow diagonal that would be invisible to physical quadrature.
    lv=logvar(yloc)
    tau=mp.sqrt(e+mp.exp(lv))
    L=l*(R/yloc)**p
    local=yloc/(32*max(1,L*p))
    if not e and tau<local/64:
        def diag(t):
            ratio=-L*mp.expm1(-p*mp.log1p(tau*t/yloc))
            u=yloc+tau*t
            return q0*mp.exp(-ratio/2-t*t*mp.exp(-ratio)/2-side*a*u-u*u/2)
        central=mp.quad(diag,[-16,-7,-2,0,2,7,16])
        left=mp.quad(f,panels(0,yloc-local,ext+(yloc/2,)))
        right=mp.quad(f,panels(yloc+local,R,ext+(yloc*1.5,yloc*2)))
        return central+left+right
    near=tuple(yloc+j*tau for j in (-16,-7,-2,0,2,7,16))
    return mp.quad(f,panels(0,R,ext+near))


def double_branch(a,y,e,p,q,l,side,k=0):
    yloc=side*y
    q0=1/math.sqrt(2*math.pi)
    lB=math.log((q*R)**2)
    def f(u):
        if u<=0:
            return 0.0
        lv=lB-l*((R/u)**p-1)
        if e:
            lv=math.log(e+ (math.exp(lv) if lv>-745 else 0.))
        dist=u-yloc
        logarg=2*math.log(abs(dist))-lv if dist else -math.inf
        if logarg>20:
            return 0.
        arg=math.exp(logarg)/2
        tot=-(k+1)*lv/2-arg-side*a*u-u*u/2
        if tot<-720:
            return 0.
        z=math.copysign(math.sqrt(2*arg),dist)
        return q0*math.exp(tot)*hermite(k,z)*side**k
    extra = [R/16,R/8,R/4,R/2,R*3/4,R*7/8]
    if 0<yloc<R:
        lv=lB-l*((R/yloc)**p-1)
        tau=math.sqrt(e+(math.exp(lv) if lv>-745 else 0))
        extra += [yloc+j*tau for j in (-16,-7,-2,0,2,7,16)]
    return quad(f,0,R,points=panels(0,R,extra)[1:-1],epsabs=4e-13,epsrel=3e-12,limit=100)[0]


def case(row,smoke):
    index,*params=row
    a,p,q,l,qq,ll=map(mp.mpf,params)
    sides = [(1,q,l)] + ([(-1,qq,ll)] if qq else [])
    js=[mp.fsum(positive_mp(a,0,0,p,side_q,side_l,sign,k) for sign,side_q,side_l in sides)
        for k in range(JETS+1)]
    jet_cross=max(abs(float(js[k])-sum(
        double_branch(float(a),0,0,float(p),float(qq),float(ll),sgn,k) for sgn,qq,ll in sides))
        /max(1,abs(float(js[k]))) for k in range(JETS+1))
    if jet_cross>2e-11:
        raise ArithmeticError((row,"jets",jet_cross))
    entries=[]
    for j in ((2,7) if smoke else DIVISIONS):
        for sign in (-1,1):
            y=mp.mpf(sign)*R/mp.mpf(2)**j
            jet=mp.fsum(z*y**k/mp.factorial(k) for k,z in enumerate(js))
            m_zero=mp.fsum(positive_mp(a,y,0,p,side_q,side_l,sgn) for sgn,side_q,side_l in sides)
            if not qq:
                m_zero+=baseline(a,y,0,-1)
            residual=m_zero-baseline(a,y,0)
            e=(y*(mp.mpf(".5") if j%2==0 else mp.mpf(".125")))**2
            # Normalize by g(a) throughout. The zero-variance plateau is exact.
            finite=[]
            for sgn, side_q, side_l in sides:
                main=positive_mp(a,y,e,p,side_q,side_l,sgn)
                check=double_branch(float(a),float(y),float(e),float(p),float(side_q),float(side_l),sgn)
                err=abs(float(main)-check)
                if err>2e-11:
                    raise ArithmeticError((row,j,sign,"finite",err))
                finite.append((main,err))
            finite_total=mp.fsum(main for main,_ in finite)
            if not qq:
                finite_total+=baseline(a,y,e,-1)
            finite_residual=finite_total-baseline(a,y,e)
            entries.append({"y":float(y),"eps":float(e),"jet_prediction":mp.nstr(jet,26),
                            "zero_eps_total":mp.nstr(m_zero,26),
                            "zero_eps_excess":mp.nstr(residual,26),
                            "zero_eps_jet_residual":mp.nstr(residual-jet,26),
                            "finite_eps_excess":mp.nstr(finite_residual,26),
                            "finite_eps_cross_max":max(err for _,err in finite)})
    center_eps=(mp.mpf(2)**-23)**2
    center=mp.fsum(positive_mp(a,0,center_eps,p,side_q,side_l,sgn) for sgn,side_q,side_l in sides)
    if not qq:
        center+=baseline(a,0,center_eps,-1)
    expected_ratio=(1+js[0])*mp.exp(-a*a/4)*mp.sqrt(2)
    return {"i":index,"a":float(a),"p":float(p),"q_plus":float(q),"ell_plus":float(l),
            "q_minus":float(qq),"ell_minus":float(ll),"jets":[mp.nstr(x,30) for x in js],
            "jet_cross_scaled_max":jet_cross,"center_limit_density_ratio_to_N_0_2":mp.nstr(expected_ratio,25),
            "center_regularized_minus_limit":mp.nstr(center-(1+js[0]),25),"rows":entries}


def worker(rows,smoke):
    mp.mp.dps=DPS
    start=utc()
    cpu=time.process_time()
    results=[case(row,smoke) for row in rows]
    return {"pid":os.getpid(),"start":start,"finish":utc(),
            "cpu_s":time.process_time()-cpu,"cases":results}


def ps(pid):
    lines=subprocess.check_output(["/bin/ps","-axo","pid=,ppid=,%cpu=,time=,state="],text=True).splitlines()
    return [r for line in lines if (r:=line.split())[1]==str(pid)]


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--cases",type=int,default=144)
    parser.add_argument("--out",default="flat-core-results.json")
    parser.add_argument("--smoke",action="store_true")
    parser.add_argument("--start-utc")
    args=parser.parse_args()
    rows=list(grid(args.cases))
    if args.smoke:
        rows=rows[:4]
    if args.start_utc:
        target=datetime.fromisoformat(args.start_utc)
        while (remain:=(target-datetime.now(timezone.utc)).total_seconds())>0:
            time.sleep(min(remain,.2))
    workers=4 if args.smoke else 12
    pid=os.getpid()
    start=utc()
    wall=time.monotonic()
    samples=[]
    print(json.dumps({"pid":pid,"start":start,"cases":len(rows),"workers":workers}),flush=True)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        batches=[[] for _ in range(workers)]
        for i,row in enumerate(rows):
            batches[(i+i//workers)%workers].append(row)
        fs=[pool.submit(worker,batch,args.smoke) for batch in batches]
        while not all(f.done() for f in fs):
            top=subprocess.check_output(["/usr/bin/top","-l","2","-s","1","-n","0"],text=True)
            sample={"utc":utc(),"whole_cpu_user_sys_idle":list(map(float,CPU_RE.findall(top)[-1])),
                    "ps":ps(pid)}
            samples.append(sample)
            print(json.dumps(sample),flush=True)
            time.sleep(.1)
        shards=[f.result() for f in fs]
    ended=utc()
    workpids={x["pid"] for x in shards}
    active=set(map(int,subprocess.check_output(["/bin/ps","-axo","pid="],text=True).split()))
    assert len(workpids)==workers and not workpids.intersection(active)
    data={"parent":pid,"start":start,"finish":ended,"wall_s":time.monotonic()-wall,
          "worker_exit_ps_verified":True,"shards":shards,"samples":samples,
          "seed":SEED,"radius":R,"dps":DPS,"python":sys.version,"platform":platform.platform(),
          "numpy":np.__version__,"scipy":scipy.__version__,"mpmath":mp.__version__}
    Path(args.out).write_text(json.dumps(data,indent=2)+"\n")
    print(json.dumps({"finish":ended,"wall_s":data["wall_s"],"workers_exited":True}),flush=True)


if __name__=="__main__":
    main()
