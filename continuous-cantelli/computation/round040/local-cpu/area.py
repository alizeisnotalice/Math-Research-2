"""Near-critical compact mass and tail area by independent adaptive quadratures."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import json, math, os, platform, re, subprocess, sys, time
from pathlib import Path
import mpmath as mp
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.stats import qmc

SEED=400118
EXPS=(7,9,12,15,18,22)
CUTS=[-32,-16,-8,-4,-2,0,2,4,8,16,32]
R=.25
PREC=47
CPU_RE=re.compile(r"CPU usage: (.*?)% user, (.*?)% sys, (.*?)% idle")


def now():
    return datetime.now(timezone.utc).isoformat(timespec="milliseconds")


def models(n):
    pts=qmc.Sobol(4,scramble=True,seed=SEED).random_base2((n-1).bit_length())[:n]
    for i,(c,b,d,j) in enumerate(pts):
        yield i,float(.6+2.4*c),float(.6*b-.3),float(.2+.4*d),float(6*j-3) if i%4 else 0.


def quadratures(profile,T,high):
    f=mp if high else math
    _,c,b,d,delta=(mp.mpf(x) for x in profile) if high else profile
    A=f.exp(1/(2*c))/f.sqrt(c)
    B1=-b*A*(c**-3+3*c**-2)/2
    K=2*A
    L=K*K/4+2*B1+delta
    q0=1/f.sqrt(2*f.pi)
    if high:
        integrate=lambda fn,cuts: mp.quad(fn,cuts)
    else:
        integrate=lambda fn,cuts: quad(fn,cuts[0],cuts[-1],points=cuts[1:-1],
                                      epsabs=8e-13,epsrel=8e-13,limit=100)[0]
    def compact(s):
        u=s/T
        return q0*f.exp(s-(c*s*s+b*s*s*u+d*s*s*u*u+u*u)/2)
    def tail(w):
        x=T+w
        perturb=(K+L/x)/x**3
        return -q0*f.exp(-w*w/2)*f.expm1(-T*T*perturb/2)
    def area(x_over_T):
        x=T*x_over_T
        p=(K+L/x)/x**3
        return -T*f.expm1(-f.log1p(p)/2-x*x*p/(2*(1+p)))
    def compact_area(s):
        u=s/T
        m=s*s*(c+b*u+d*u*u)
        if abs(s)<(mp.mpf("1e-19") if high else 1e-7):
            return q0*(1+1.5*s)
        sig=f.sqrt(m)
        z1=(sig-s/sig)/f.sqrt(2)
        z2=z1+sig/f.sqrt(2)
        if z2<0:
            diff=(f.erfc(-z2)-f.erfc(-z1))/2
        else:
            diff=(f.erfc(z1)-f.erfc(z2))/2
        return q0*f.sqrt(2*f.pi)*f.exp(s*s/(2*m)-u*u/2)*diff/sig

    C=integrate(compact,CUTS)
    D=T*integrate(tail,CUTS)
    X=integrate(area,[1,1.25,1.5,1.75,2])
    Y=integrate(compact_area,CUTS)
    return C,D,X,Y,A,B1,K,L,delta


def work(batch,exps):
    mp.mp.dps=PREC
    start=now(); cpu=time.process_time(); models_out=[]
    for profile in batch:
        rows=[]
        for j in exps:
            T=2**j
            hi=quadratures(profile,mp.mpf(T),True)
            lo=quadratures(profile,float(T),False)
            C,D,X,Y,A,A1,K,L,delta=hi
            cross=max(abs(float(x)-y) for x,y in zip(hi[:4],lo[:4]))
            if cross>5e-11:
                raise ArithmeticError((profile,j,cross))
            rows.append({"exp":j,"point":mp.nstr(T*(D-C),25),"point_limit":mp.nstr(delta/2,25),
                         "area":mp.nstr(T*(X-Y),25),"area_limit":mp.nstr(delta/4,25),
                         "area_vs_Alog_two":mp.nstr(T*(X-A*mp.log(2)),25),
                         "comp_area_vs_Alog_two":mp.nstr(T*(Y-A*mp.log(2)),25),
                         "comp_area_limit":mp.nstr(A1/2,25),
                         "scaled_point_difference_error":float(abs(T*(D-C)-delta/2)),
                         "scaled_area_difference_error":float(abs(T*(X-Y)-delta/4)),
                         "crosscheck":cross})
        models_out.append({"model":profile,"K":mp.nstr(K,25),"L":mp.nstr(L,25),"rows":rows})
    return {"pid":os.getpid(),"start":start,"finish":now(),"cpu_seconds":time.process_time()-cpu,"models":models_out}


def main():
    args=argparse.ArgumentParser()
    args.add_argument("--cases",type=int,default=256)
    args.add_argument("--out",default="area-results.json")
    args.add_argument("--start-utc")
    args.add_argument("--smoke",action="store_true")
    a=args.parse_args()
    data=list(models(a.cases))
    n=2 if a.smoke else 12
    exps=(7,22) if a.smoke else EXPS
    if a.start_utc:
        launch=datetime.fromisoformat(a.start_utc)
        while (remain:=(launch-datetime.now(timezone.utc)).total_seconds())>0:
            time.sleep(min(remain,.2))
    start=now(); pid=os.getpid(); wall=time.monotonic(); samples=[]
    print(json.dumps({"parent":pid,"start":start,"workers":n,"profiles":len(data)}),flush=True)
    with ProcessPoolExecutor(n) as pool:
        tasks=[pool.submit(work,data[i::n],exps) for i in range(n)]
        while not all(f.done() for f in tasks):
            stat=subprocess.check_output(["/usr/bin/top","-l","2","-s","1","-n","0"],text=True)
            pids=subprocess.check_output(["/bin/ps","-axo","pid=,ppid=,%cpu=,time=,state="],text=True)
            sample={"utc":now(),"machine_cpu_user_sys_idle":list(map(float,CPU_RE.findall(stat)[-1])),
              "workers":[r for line in pids.splitlines() if (r:=line.split())[1]==str(pid)]}
            samples.append(sample); print(json.dumps(sample),flush=True); time.sleep(.1)
        vals=[task.result() for task in tasks]
    stop=now(); pids={r["pid"] for r in vals}
    active={int(v) for v in subprocess.check_output(["/bin/ps","-axo","pid="],text=True).split()}
    assert len(pids)==n and not pids.intersection(active)
    result={"start":start,"finish":stop,"wall_s":time.monotonic()-wall,"parent":pid,
        "worker_exit_verified":True,"shards":vals,"samples":samples,
        "dps":PREC,"seed":SEED,"python":sys.version,"platform":platform.platform(),
        "mpmath":mp.__version__,"numpy":np.__version__,"scipy":scipy.__version__}
    Path(a.out).write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"finish":stop,"wall_s":result["wall_s"],"workers_exited":True}),flush=True)


if __name__=="__main__":
    main()
