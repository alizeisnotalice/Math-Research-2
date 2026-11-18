"""First-hitting heat potentials for slowly decaying infinite-order flatness."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import json, math, os, re, subprocess, sys, time
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.special import erfcx
from scipy.stats import qmc

SEED=360118
LEVELS=(2,5,10,20,45,90,180,300)
DPS=48
RADIUS=.25
KAPPA=4


def now():
    return datetime.now(timezone.utc).isoformat(timespec="milliseconds")


def config(n):
    S=qmc.Sobol(5,scramble=True,seed=SEED).random_base2((n-1).bit_length())[:n]
    for i,(p,gap,a,flat,minus) in enumerate(S):
        yield i,float(1.15+1.3*p),float(.75+3.75*gap),float(6*a-2),float(100+120*flat),float(.25+.75*minus)


def geometry(profile,j,m):
    _,p,A0,a,q0,minus=profile
    L0=m(KAPPA)
    d=(A0+2)/(p*L0**(p-1))
    r0=m(RADIUS)
    v0=r0*r0/q0
    Lr=(L0**p+(j-1)*m(2).log()/d)**(1/p) if m is not float else (L0**p+(j-1)*math.log(2)/d)**(1/p)
    Ly=(L0**p+j*m(2).log()/d)**(1/p) if m is not float else (L0**p+j*math.log(2)/d)**(1/p)
    return p,d,a,minus,Lr,Ly,v0*2**(-(j-1))


def mp_parameters(profile,j):
    _,p,A0,a,q0,minus=(mp.mpf(s) for s in profile)
    L0=mp.mpf(KAPPA)
    d=(A0+2)/(p*L0**(p-1))
    r0=mp.mpf(RADIUS)
    Lr=(L0**p+(j-1)*mp.log(2)/d)**(1/p)
    Ly=(L0**p+j*mp.log(2)/d)**(1/p)
    return p,d,a,minus,Lr,Ly,r0*r0/q0*mp.mpf(2)**(-(j-1))


def W_mp(q):
    if q<10**11:
        return q*(1-mp.sqrt(mp.pi*q/2)*mp.exp(q/2)*mp.erfc(mp.sqrt(q/2)))
    return 1-3/q+15/q**2-105/q**3+945/q**4


def W_float(q):
    if q<8000:
        return q*(1-math.sqrt(math.pi*q/2)*float(erfcx(math.sqrt(q/2))))
    t=1/q
    return 1+t*(-3+t*(15+t*(-105+t*(945-10395*t))))


def evaluate(profile,j,high):
    f=mp if high else math
    p,d,a,minus,Lr,Ly,vR=mp_parameters(profile,j) if high else geometry(profile,j,float)
    R=RADIUS*f.exp(KAPPA-Lr)
    y=RADIUS*f.exp(KAPPA-Ly)
    h=f.sqrt(vR/2)
    qr=R*R/vR
    A=d*p*Lr**(p-1)-2
    qgauss=1/f.sqrt(2*f.pi)
    W=W_mp if high else W_float
    if high:
        integrate=lambda g,xs: mp.quad(g,xs)
    else:
        integrate=lambda g,xs: quad(g,xs[0],xs[-1],points=xs[1:-1],
                                     epsabs=4e-13,epsrel=3e-12,limit=80)[0]

    def center(t):
        dx=2*t/(A*qr)
        delta=d*Lr**p*f.expm1(p*f.log1p(dx/Lr))
        qshift=delta-2*dx
        u=R*f.exp(-dx)
        sdiff=-a*(u-R)-(u*u-R*R)/2
        return W(qr*f.exp(qshift))*f.exp(-qr*f.expm1(qshift)/2-qshift-delta/2-dx+sdiff)

    central=integrate(center,[0,1,4,10,20,35,50])

    def hit(t):
        ru=1+h*t/y
        if not ru:
            return 0
        dl=-f.log1p(h*t/y) if ru>0 else -f.log(-ru)
        logva=-d*Ly**p*f.expm1(p*f.log1p(dl/Ly))
        if ru<0:
            logva+=f.log(minus)
        if logva<-32:
            return 0
        var=f.exp(logva)
        z=abs(t)/f.sqrt(var)
        if z>30:
            return 0
        heat=f.sqrt(var)*(qgauss*f.exp(-z*z/2)-z*f.erfc(z/f.sqrt(2))/2)
        ht=h*t
        return 2*heat*f.exp(-a*ht-y*ht-ht*ht/2)

    lo=max(-24,(-R-y)/h)
    hi=min(24,(R-y)/h)
    cuts=[lo]+[x for x in (-12,-6,-2,0,2,6,12) if lo<x<hi]+[hi]
    mass=integrate(hit,cuts)
    log_pref=f.log(qgauss)+f.log(R)+f.log(vR)/2+f.log(2)-f.log(A)-2*f.log(qr)-a*R-R*R/2
    log_pos=log_pref-qr/2+f.log(central)
    return central,mass,log_pos,R,y,qr,A


def worker(batch):
    mp.mp.dps=DPS
    begin=now(); cpu=time.process_time()
    results=[]
    for profile in batch:
        rows=[]
        for j in LEVELS:
            hp=evaluate(profile,j,True)
            dp=evaluate(profile,j,False)
            error=max(abs(float(hp[i])-dp[i]) for i in (0,1))
            if error>5e-11:
                raise ArithmeticError((profile,j,error))
            C,M,L,R,y,q,A=hp
            rows.append({"level":j,"outer_radius":mp.nstr(R,22),"hitting_radius":mp.nstr(y,22),
                         "relative_radius_gap":mp.nstr(1-y/R,22),
                         "positive_side_center_log_normalized":mp.nstr(L,25),
                         "center_log_divide_q":mp.nstr(-L/q,22),
                         "q":mp.nstr(q,22),"q_slope":mp.nstr(A,22),
                         "normalized_endpoint_integral":mp.nstr(C,25),
                         "normalized_hitting_mass":mp.nstr(M,25),
                         "analytic_jet_count_budget":mp.nstr(q/(8*mp.log(4)),22),
                         "crosserror":error})
        results.append({"profile":profile,"rows":rows})
    return {"pid":os.getpid(),"start":begin,"finish":now(),
            "cpu_seconds":time.process_time()-cpu,"profiles":results}


CPU_RE=re.compile(r"CPU usage: (.*?)% user, (.*?)% sys, (.*?)% idle")


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--cases",type=int,default=320)
    parser.add_argument("--out",default="heat-results.json")
    parser.add_argument("--workers",type=int,default=12)
    args=parser.parse_args()
    profiles=list(config(args.cases)); parent=os.getpid()
    start=now(); begin=time.monotonic(); samples=[]
    print(json.dumps({"pid":parent,"start":start,"workers":args.workers}),flush=True)
    with ProcessPoolExecutor(args.workers) as pool:
        futures=[pool.submit(worker,profiles[j::args.workers]) for j in range(args.workers)]
        while not all(f.done() for f in futures):
            supply=subprocess.check_output(["/usr/bin/top","-l","2","-s","1","-n","0"],text=True)
            processes=subprocess.check_output(["/bin/ps","-axo","pid=,ppid=,%cpu=,time=,state="],text=True)
            sample={"utc":now(),"machine_cpu_user_sys_idle":list(map(float,CPU_RE.findall(supply)[-1])),
                    "workers":[r for line in processes.splitlines() if (r:=line.split())[1]==str(parent)]}
            samples.append(sample); print(json.dumps(sample),flush=True); time.sleep(.1)
        shards=[f.result() for f in futures]
    finish=now()
    pids={s["pid"] for s in shards}
    exists={int(v) for v in subprocess.check_output(["/bin/ps","-axo","pid="],text=True).split()}
    assert len(pids)==args.workers and not pids.intersection(exists)
    doc={"start":start,"finish":finish,"wall_s":time.monotonic()-begin,"parent":parent,
         "qualified_worker_exit":True,"samples":samples,"shards":shards,"versions":{
          "python":sys.version,"numpy":np.__version__,"scipy":scipy.__version__,"mpmath":mp.__version__},
         "precision":DPS,"seed":SEED}
    Path(args.out).write_text(json.dumps(doc,indent=2)+"\n")
    print(json.dumps({"finish":finish,"wall_s":doc["wall_s"],"workers_exited":True}),flush=True)


if __name__=="__main__":
    main()
