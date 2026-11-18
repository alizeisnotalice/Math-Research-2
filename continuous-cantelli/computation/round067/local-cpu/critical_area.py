"""Near-axis heat-potential remainder and window area for log-flat models."""
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import argparse, json, math, os, platform, re, subprocess, sys, time
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy import __version__ as scipy_version
from scipy.integrate import quad
from scipy.stats import qmc

LEVELS=(12,24,40,72,120)
J0,J1=.6,1.4
EDGE=2.
LOG_CORE=4.
DPS=53
SEED=631118
TICKS="CPU usage: (.*?)% user, (.*?)% sys, (.*?)% idle"


def utc():
    return datetime.now(timezone.utc).isoformat(timespec="milliseconds")


def profiles(n):
    data=qmc.Sobol(5,scramble=True,seed=SEED).random_base2((n-1).bit_length())[:n]
    for i,(c,offset,window,p,a) in enumerate(data):
        k=(.25,.5,1.,2.)[i%4]
        yield i,k,(i//4)%2==1,float(.22+.75*c),float(.35+.45*offset),float(2*window-1),float(.4+.8*p),float(6*a-3)


def panel(a,b,interior):
    return [a]+sorted({t for t in interior if a<t<b})+[b]


def scaled_ports(c,k,phase,narrow,L,f):
    bands={1:[],-1:[]}
    for L0 in LEVELS:
        rho=f.sqrt(c)/L0**(k/2)
        w=min(.20, .8/L0**narrow)
        for sign,center in ((1,J0*f.exp(phase*rho)),(1,J1*f.exp(-phase*rho/2)),
                            (1,2.25*f.exp(phase*.09)),(-1,2.65*f.exp(phase*.12))):
            # r0/r = exp(L-L0). Hold the five bands fixed in the physical model.
            bands[sign].append((L-L0+f.log(center),w))
    return bands


def calculate(profile,L,high,nodes,weights):
    _,k,sparse,c,y,phase,narrow,a=profile
    f=mp if high else math
    if high:
        k,c,y,phase,narrow,a,LL=(mp.mpf(v) for v in (k,c,y,phase,narrow,a,L))
    else:
        LL=float(L)
    rho=f.exp(-LL)
    q0=1/f.sqrt(2*f.pi)
    root2=f.sqrt(2)
    bands=scaled_ports(c,k,phase,narrow,LL,f) if sparse else {1:[],-1:[]}
    if high:
        def integrate(fn,cuts):
            return f.quad(fn,cuts)
    else:
        def integrate(fn,cuts):
            return quad(fn,cuts[0],cuts[-1],points=cuts[1:-1],epsabs=8e-13,epsrel=3e-11,limit=100)[0]

    def q(w,sign):
        func=c/(LL-w)**k
        if not sparse:
            return func
        shape=max(0, *(1-abs((w-w0)/width) for w0,width in bands[sign]))
        return func*shape**2

    def envelope(w,sign):
        Q=q(w,sign)
        if not Q:
            return 0
        d=y*f.exp(-w)
        x=sign*rho*f.exp(w)
        values=[]
        for t,b in zip(nodes,weights):
            v=1-sign*t*d
            power=-v*v/(2*Q)-f.log(Q)*1.5
            if power < -650:
                continue
            values.append(b*(1-t)**2*v*f.exp(power))
        return y**3*q0*f.exp(-w-a*x-x*x/2)*sum(values)

    W0=f.log(EDGE)
    Wend=LL-LOG_CORE
    rem=[]
    for sign in (1,-1):
        if sparse:
            pieces=[]
            for w,half in bands[sign]:
                lo,hi=max(W0,w-half),min(Wend,w+half)
                if lo<hi:
                    pieces.append(integrate(lambda t:envelope(t,sign),panel(lo,hi,(w,))))
            rem.append(sum(pieces))
        else:
            rem.append(integrate(lambda w:envelope(w,sign),panel(W0,Wend,(1,2,4,8,16,30,50,80))))

    def tail(t):
        if t>35:
            return 0.
        return (1+t*t)*f.erfc(t/root2)/2-t*q0*f.exp(-t*t/2)

    def window(aa,bb,sigma):
        A,B=aa/sigma,bb/sigma
        if A>0:
            return tail(A)-tail(B)
        if B<0:
            return tail(-B)-tail(-A)
        return 1-tail(-A)-tail(B)

    def density(t,correction):
        if not t:
            return 0
        sig=1 if t>0 else -1
        Q=q(f.log(abs(t)),sig)*t*t
        if not Q:
            return 0
        T=rho*t
        g=f.exp(-a*T-T*T/2)
        if not correction:
            return Q*g
        return Q*g*(window(J0-t,J1-t,f.sqrt(Q))-(J0<t<J1))

    extras=[-1,0,.25,.45,J0,.8,1.,1.2,J1,1.6]
    if sparse:
        for sig in (-1,1):
            for w,width in bands[sig]:
                for qv in (w-width,w,w+width):
                    if qv<f.log(EDGE):
                        extras.append(sig*f.exp(qv))
    sharp=integrate(lambda t:density(t,False),panel(J0,J1,extras))
    area_delta=integrate(lambda t:density(t,True),panel(-EDGE,EDGE,extras))
    return rem[0],rem[1],sharp,area_delta


def worker(batch,smoke):
    mp.mp.dps=DPS
    cpu=time.process_time(); start=utc()
    # Build independently chosen quadrature orders for the inner remainder integral.
    m_nodes,m_weights=mp.gauss_quadrature(7,"legendre")
    D_nodes,D_weights=np.polynomial.legendre.leggauss(16)
    M=tuple((x+1)/2 for x in m_nodes),tuple(x/2 for x in m_weights)
    D=tuple(float((x+1)/2) for x in D_nodes),tuple(float(x/2) for x in D_weights)
    cases=[]
    for profile in batch:
        rows=[]
        for L in ((12,120) if smoke else LEVELS):
            hi=calculate(profile,L,True,*M)
            lo=calculate(profile,L,False,*D)
            errors=[abs(float(x)-y) for x,y in zip(hi,lo)]
            if max(errors)>2e-9*max(1,*map(abs,lo)):
                raise ArithmeticError((profile,L,errors,list(map(float,hi)),lo))
            plus,minus,mass,de=hi
            natural=mp.mpf(profile[3])*L**(-mp.mpf(profile[1]))
            rows.append({"L":L,"positive_taylor":mp.nstr(plus,25),"negative_taylor_magnitude":mp.nstr(minus,25),
                         "absolute_remainder_after_natural_scaling":mp.nstr((plus+minus)/natural,25),
                         "area_sharp":mp.nstr(mass,25),"area_signed_error":mp.nstr(de,25),
                         "relative_area_error":mp.nstr(de/mass,25) if mass>0 else None,
                         "area_error_by_natural_scale":mp.nstr(de/natural,25),
                         "max_cross_error":max(errors)})
        cases.append({"profile":profile,"rows":rows})
    return {"pid":os.getpid(),"start":start,"finish":utc(),"cpu_seconds":time.process_time()-cpu,"cases":cases}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--cases",type=int,default=88)
    p.add_argument("--out",default="critical-area-results.json")
    p.add_argument("--smoke",action="store_true")
    args=p.parse_args()
    cfg=list(profiles(args.cases)); workers=2 if args.smoke else 12
    parent=os.getpid(); begin=utc(); clock=time.monotonic(); samples=[]
    print(json.dumps({"pid":parent,"start":begin,"workers":workers,"cases":len(cfg)}),flush=True)
    with ProcessPoolExecutor(workers) as pool:
        # Rotate shard assignments to balance smooth/sparse and k values.
        batches=[[] for _ in range(workers)]
        for i,row in enumerate(cfg):
            batches[(i+i//workers)%workers].append(row)
        fs=[pool.submit(worker,rows,args.smoke) for rows in batches]
        while not all(f.done() for f in fs):
            ts=subprocess.check_output(["/usr/bin/top","-l","2","-s","1","-n","0"],text=True)
            ps=subprocess.check_output(["/bin/ps","-axo","pid=,ppid=,%cpu=,time=,state="],text=True)
            s={"utc":utc(),"cpu_pct_user_sys_idle":list(map(float,re.findall(TICKS,ts)[-1])),
                "ps":[a for x in ps.splitlines() if (a:=x.split())[1]==str(parent)]}
            print(json.dumps(s),flush=True); samples.append(s); time.sleep(.1)
        outputs=[f.result() for f in fs]
    finish=utc()
    IDs={o["pid"] for o in outputs}
    active={int(v) for v in subprocess.check_output(["/bin/ps","-axo","pid="],text=True).split()}
    assert len(IDs)==workers and not IDs.intersection(active)
    out={"parent":parent,"start":begin,"finish":finish,"wall_seconds":time.monotonic()-clock,
         "shards":outputs,"samples":samples,"exit_verified":True,"dps":DPS,"seed":SEED,
         "python":sys.version,"platform":platform.platform(),"scipy":scipy_version,
         "mpmath":mp.__version__,"numpy":np.__version__}
    Path(args.out).write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps({"finish":finish,"wall_seconds":out["wall_seconds"],"exit_verified":True}),flush=True)


if __name__=="__main__":
    main()
