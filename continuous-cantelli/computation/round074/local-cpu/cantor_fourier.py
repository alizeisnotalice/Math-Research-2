"""Complementary-gap Fourier coefficients of distance-powered Cantor variance."""
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import argparse, cmath, json, math, os, platform, re, subprocess, sys, time
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy import __version__ as scipy_version
from scipy.integrate import quad

DPS=46
TOL=2e-18
FAMILIES=((3,(0,2)),(4,(0,3)),(5,(0,2,4)),(5,(0,1,4)),
          (7,(0,2,4,6)),(8,(0,1,3,7)))
ALPHAS=(.5,.8,1.2,1.65,2.,2.35,2.7,3.1)
TICKS="CPU usage: (.*?)% user, (.*?)% sys, (.*?)% idle"


def utc():
    return datetime.now(timezone.utc).isoformat(timespec="milliseconds")


def profiles(n):
    for i in range(n):
        b,D=FAMILIES[i%len(FAMILIES)]
        a=ALPHAS[(i//len(FAMILIES))%len(ALPHAS)]
        yield i,b,D,a,.09+.88*((13*i+5)%97)/96


def gaps(base,digits):
    return tuple((i+1,j) for i,j in zip(digits,digits[1:]) if j>i+1)


def local_integrals(theta,mu,a,high,rule):
    f=mp if high else math
    if mu<=30:
        p=math.ceil(4/float(a))
        parts=max(1,math.ceil(float(theta)/15))

        def integrand(t,cosine):
            x=t**p
            osc=f.cos(theta*(1-x)) if cosine else 1
            return p*t**(p-1)*osc*f.expm1(-mu*x**a)
        return tuple(integrate(lambda t:integrand(t,j),parts,high,rule) for j in (True,False))

    # Subtract the unattenuated cosine integral analytically, locating the
    # sharply damped endpoints with u = mu*x^a.
    cap=min(mu,100)
    exp=math.ceil(4*float(a))
    p=exp/a
    scale=(cap/mu)**(1/a)
    parts=max(1,math.ceil(float(theta*scale)/12))

    def integrand(t,cosine):
        u=cap*t**exp
        x=scale*t**p
        osc=f.cos(theta*(1-x)) if cosine else 1
        return scale*p*t**(p-1)*osc*f.exp(-u)

    return (integrate(lambda t:integrand(t,True),parts,high,rule)-f.sin(theta)/theta,
            integrate(lambda t:integrand(t,False),parts,high,rule)-1)


def integrate(fn,parts,high,rule):
    if high:
        nodes,weights=rule
        step=mp.mpf(1)/parts
        half=step/2
        return half*mp.fsum(w*fn((j+mp.mpf(".5"))*step+half*x)
                            for j in range(parts) for x,w in zip(nodes,weights))
    cuts=[j/parts for j in range(1,parts)]
    return quad(fn,0,1,points=cuts,epsabs=3e-13,epsrel=3e-12,
                limit=max(90,8*parts))[0]


def calc(base,digits,alpha,c,n,high,rule=None,fixed_depth=None):
    f=mp if high else math
    a=f.mpf(alpha) if high else alpha
    c=f.mpf(c) if high else c
    q=len(digits)/base
    ratio=len(digits)/base**(float(a)+1)
    base_gaps=gaps(base,digits)
    damp=2*f.pi*f.pi*n*n*c
    prefactor=float(damp)*sum(((j-i)/base)**(1+float(a))
                              for i,j in base_gaps)/((1+float(a))*2**float(a))
    depth=max(5,math.ceil(math.log(TOL*(1-ratio)/prefactor)/math.log(ratio)))
    if fixed_depth is not None:
        depth=fixed_depth
    bound=prefactor*ratio**depth/(1-ratio)
    scale=f.mpf(1) if high else 1.
    E=mp.mpc(1) if high else 1+0j
    ex=mp.exp if high else cmath.exp
    pi=f.pi
    total=mp.mpc(0) if high else 0j
    mass=f.mpf(1) if high else 1.
    denom=base
    for j in range(depth):
        for lo,hi in base_gaps:
            width=(f.mpf(hi-lo) if high else hi-lo)/denom
            theta=pi*n*width
            mu=damp*(width/2)**a
            val,m=local_integrals(theta,mu,a,high,rule)
            weight=(f.mpf(hi-lo) if high else hi-lo)/base*scale
            center=ex(1j*pi*((n*(lo+hi))%(2*denom))/denom)
            total+=E*center*weight*val
            mass+=weight*m
        E*=sum(ex(2j*pi*((n*d)%denom)/denom) for d in digits)/len(digits)
        scale*=f.mpf(len(digits))/base if high else q
        denom*=base
    return total,mass,depth,bound


def cantor_phase(base,digits):
    phase=mp.mpc(1)
    denom=base
    for j in range(95):
        phase*=sum(mp.exp(2j*mp.pi*d/denom) for d in digits)/len(digits)
        denom*=base
    return phase


def prefix_check(base,digits,n):
    V=[0]
    prod=mp.mpc(1)
    errs=[]
    denom=1
    for depth in range(1,5):
        V=[base*x+a for x in V for a in digits]
        denom*=base
        rhs=mp.fsum(mp.exp(2j*mp.pi*((n*x)%denom)/denom) for x in V)/len(V)
        prod*=mp.fsum(mp.exp(2j*mp.pi*((n*d)%denom)/denom) for d in digits)/len(digits)
        errs.append(abs(rhs-prod))
    return max(errs)


def fmt(x):
    return [mp.nstr(mp.re(x),28),mp.nstr(mp.im(x),28)]


def worker(data,smoke,deadline):
    mp.mp.dps=DPS
    started=utc(); cpu=time.process_time()
    X,W=mp.gauss_quadrature(38,"legendre")
    rule=tuple(X),tuple(W)
    fine=None
    tests=[]
    for i,b,digits,alpha,c in data:
        if deadline is not None and time.monotonic()>=deadline and tests:
            break
        exps=(2,5 if alpha>2 else 6) if smoke else (2,4,5 if alpha>2 else 6)
        ns=[b**s for s in exps]
        if not smoke:
            ns.append(b**4+1)
        phase=cantor_phase(b,digits)
        rows=[]
        for n in ns:
            h,m,depth,tail=calc(b,digits,alpha,c,n,True,rule)
            d,dm,_,_=calc(b,digits,alpha,c,n,False)
            cross=max(abs(complex(h)-d),abs(float(m)-dm))
            if cross>1e-9:
                if fine is None:
                    xx,ww=mp.gauss_quadrature(68,"legendre")
                    fine=tuple(xx),tuple(ww)
                h,m,depth,tail=calc(b,digits,alpha,c,n,True,fine)
                cross=max(abs(complex(h)-d),abs(float(m)-dm))
            rows.append({"n":n,"resonant":n in [b**s for s in exps],"coefficient":fmt(h),
                         "damped_mass":mp.nstr(m,28),"coefficient_divided_by_mass":fmt(h/m),
                         "phase_ratio_error":mp.nstr(abs(h/m-phase),25),
                         "cross_difference":cross,"cross_pass":cross<=1e-9,
                         "prefix_four_levels_error":mp.nstr(prefix_check(b,digits,n),14),
                         "omitted_level_start":depth,"tail_bound":tail})
        tests.append({"profile":i,"base":b,"digits":digits,"alpha":alpha,"c":c,
                      "cantor_phase_limit":fmt(phase),"rows":rows})
    return {"pid":os.getpid(),"start":started,"finish":utc(),
            "cpu_seconds":time.process_time()-cpu,"tests":tests}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--cases",type=int,default=48); p.add_argument("--out",default="cantor-results.json")
    p.add_argument("--smoke",action="store_true")
    p.add_argument("--time-budget",type=float)
    a=p.parse_args()
    cfg=list(profiles(a.cases)); workers=2 if a.smoke else 12
    parent=os.getpid(); start=utc(); wall=time.monotonic(); samples=[]
    print(json.dumps({"pid":parent,"start":start,"workers":workers,"cases":len(cfg)}),flush=True)
    with ProcessPoolExecutor(workers) as pool:
        bags=[[] for _ in range(workers)]
        for j,v in enumerate(cfg):
            bags[(j+j//workers)%workers].append(v)
        deadline=wall+a.time_budget if a.time_budget else None
        fs=[pool.submit(worker,bag,a.smoke,deadline) for bag in bags]
        while not all(f.done() for f in fs):
            t=subprocess.check_output(["/usr/bin/top","-l","2","-s","1","-n","0"],text=True)
            r=subprocess.check_output(["/bin/ps","-axo","pid=,ppid=,%cpu=,time=,state="],text=True)
            sample={"utc":utc(),"machine_cpu_pct":list(map(float,re.findall(TICKS,t)[-1])),
                    "ps":[parts for line in r.splitlines() if (parts:=line.split())[1]==str(parent)]}
            samples.append(sample); print(json.dumps(sample),flush=True); time.sleep(.1)
        out=[f.result() for f in fs]
    finish=utc()
    pids={s["pid"] for s in out}
    active={int(x) for x in subprocess.check_output(["/bin/ps","-axo","pid="],text=True).split()}
    assert len(pids)==workers and not pids.intersection(active)
    doc={"start":start,"finish":finish,"parent":parent,"wall_s":time.monotonic()-wall,"worker_exit_verified":True,
         "samples":samples,"shards":out,"time_budget_s":a.time_budget,"mp_dps":DPS,"python":sys.version,"platform":platform.platform(),
         "mpmath":mp.__version__,"scipy":scipy_version,"numpy":np.__version__}
    Path(a.out).write_text(json.dumps(doc,indent=2)+"\n")
    print(json.dumps({"finish":finish,"wall_s":doc["wall_s"],"worker_exit_verified":True}),flush=True)


if __name__=="__main__":
    main()
