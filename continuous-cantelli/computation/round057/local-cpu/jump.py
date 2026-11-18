"""Diagonal Gaussian density jet for asymmetric one-sided contact."""
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import argparse, json, math, os, platform, re, subprocess, sys, time
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy import __version__ as scipy_version
from scipy.integrate import quad
from scipy.stats import qmc

R=.25
DPS=87
SEED=521118
ORDERS=((3,3),(3,4),(4,4),(4,5),(5,6),(6,6),(6,8),(8,8),(10,12))
PANELS=(0,1,3,6,10,15,20)
CPU_RE=re.compile(r"CPU usage: (.*?)% user, (.*?)% sys, (.*?)% idle")


def utc():
    return datetime.now(timezone.utc).isoformat(timespec="milliseconds")


def grid(n):
    rows=qmc.Sobol(7,scramble=True,seed=SEED).random_base2((n-1).bit_length())[:n]
    for j,(a,c1,c2,b1,b2,d1,d2) in enumerate(rows):
        yield j,float(6*a-3),ORDERS[j%len(ORDERS)],(
          [float(.25+2.25*c1),float(1.2*b1-.6),float(.1+.25*d1)],
          [float(.25+2.25*c2),float(1.2*b2-.6),float(.1+.25*d2)])


def high(n,a,c,b,d,z,sig):
    a,c,b,d=(mp.mpf(q) for q in (a,c,b,d))
    z=mp.mpf(z)
    H=1+b*z+d*z*z
    e=mp.sqrt(c*z**(n-2)*H)
    q0=1/mp.sqrt(2*mp.pi)

    def term(t):
        u=z*e*t
        H_ratio=((b+2*d*z)*u+d*u*u)/H
        # Only local perturbation of the log-variance enters.
        lv=n*mp.log1p(e*t)+mp.log1p(H_ratio)
        delta=-sig*(a+sig*z)*u-u*u/2-lv/2-t*t*mp.expm1(-lv)/2
        return mp.expm1(delta)

    raw=mp.quad(lambda t:q0*mp.exp(-t*t/2)*(term(t)+term(-t)),PANELS)
    return mp.exp(-sig*a*z-z*z/2)*raw/z**(n-2)


def fast(n,a,c,b,d,z,sig):
    H=1+b*z+d*z*z
    e=math.sqrt(c*z**(n-2)*H)
    q0=1/math.sqrt(2*math.pi)

    def paired(t):
        et=e*t
        u=z*et
        perturb=(b+2*d*z)*u/H
        even=d*u*u/H
        lv_mean=(n*math.log1p(-et*et)+math.log1p(2*even+even*even-perturb*perturb))/2
        lv_half=n*math.atanh(et)+math.atanh(perturb/(1+even))
        exp_lv=math.exp(-lv_mean)
        mean=-u*u/2-lv_mean/2-t*t*(math.expm1(-lv_mean)+2*exp_lv*math.sinh(lv_half/2)**2)/2
        half=-sig*(a+sig*z)*u-lv_half/2+t*t*exp_lv*math.sinh(lv_half)/2
        both=2*(math.expm1(mean)+2*math.exp(mean)*math.sinh(half/2)**2)
        return q0*math.exp(-t*t/2)*both/z**(n-2)

    integral=quad(paired,PANELS[0],PANELS[-1],points=PANELS[1:-1],
                  epsabs=1e-10,epsrel=3e-12,limit=100)[0]
    return math.exp(-sig*a*z-z*z/2)*integral


def leading(n,a,c,b,d,sig):
    C0=mp.mpf(c)*n*(n-1)/2
    C1=mp.mpf(c)*n*(n+1)*(mp.mpf(b)-sig*mp.mpf(a))/2
    if n==3:
        C1+=45*mp.mpf(c)**2
    return C0,C1


def run(rows,smoke):
    mp.mp.dps=DPS
    start=utc(); cpu=time.process_time()
    result=[]
    for idx,a,orders,shape in rows:
        n0=min(orders)
        exps=(11,13,15,17,19,21) if n0==3 else ((6,8,10,13,16,18) if n0<8 else (4,6,8,10,12,14))
        if smoke:
            exps=exps[0],exps[-1]
        info=[]
        C=[leading(n,a,*data,sig) for n,data,sig in zip(orders,shape,(1,-1))]
        fac=mp.factorial(n0-2)
        target=fac*sum(mp.mpf(sig)**(n0+1)*C[j][0] for j,(n,sig) in enumerate(zip(orders,(1,-1))) if n==n0)
        for ex in exps:
            z=R*2**(-ex)
            halves=[]
            for n, data, sig in zip(orders,shape,(1,-1)):
                h=high(n,a,*data,z,sig)
                f=fast(n,a,*data,z,sig)
                err=abs(float(h)-f)
                if err>3e-10*max(1,abs(f)):
                    raise ArithmeticError((idx,n,ex,str(h),f,err))
                halves.append((h,err))
            J=fac*sum(mp.mpf(sig)**(n0+1)*halves[k][0]*mp.mpf(z)**(n-n0)
                      for k,(n,sig) in enumerate(zip(orders,(1,-1))))
            info.append({"exp":ex,"z":z,"observed_jump":mp.nstr(J,29),
                         "jump_error":mp.nstr(J-target,25),
                         "sides":[{"sig":sig,"order":n,"normalized":mp.nstr(h,33),
                                   "expected":mp.nstr(C[j][0],26),"first_slope":mp.nstr(C[j][1],26),
                                   "scaled_slope_error":mp.nstr((h-C[j][0])/z-C[j][1],22),
                                   "cross_error":err}
                                  for j,((h,err),n,sig) in enumerate(zip(halves,orders,(1,-1)))]})
        result.append({"index":idx,"center":a,"orders":orders,"coefficients":shape,
                       "predicted_jump_divided_by_source":mp.nstr(target,28),"rows":info})
    return {"pid":os.getpid(),"start":start,"finish":utc(),
            "cpu_seconds":time.process_time()-cpu,"cases":result}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--cases",type=int,default=144)
    p.add_argument("--out",default="jump-results.json")
    p.add_argument("--smoke",action="store_true")
    a=p.parse_args()
    data=list(grid(a.cases)); workers=2 if a.smoke else 12
    pid=os.getpid(); start=utc(); wall=time.monotonic(); samples=[]
    print(json.dumps({"parent":pid,"start":start,"workers":workers,"cases":len(data)}),flush=True)
    with ProcessPoolExecutor(workers) as pool:
        fs=[pool.submit(run,data[i::workers],a.smoke) for i in range(workers)]
        while not all(f.done() for f in fs):
            top=subprocess.check_output(["/usr/bin/top","-l","2","-s","1","-n","0"],text=True)
            table=subprocess.check_output(["/bin/ps","-axo","pid=,ppid=,%cpu=,time=,state="],text=True)
            row={"utc":utc(),"whole_cpu_user_sys_idle":list(map(float,CPU_RE.findall(top)[-1])),
                 "ps":[r for ln in table.splitlines() if (r:=ln.split())[1]==str(pid)]}
            samples.append(row); print(json.dumps(row),flush=True); time.sleep(.1)
        out=[f.result() for f in fs]
    finish=utc()
    pids={x["pid"] for x in out}
    active={int(v) for v in subprocess.check_output(["/bin/ps","-axo","pid="],text=True).split()}
    assert len(pids)==workers and not pids.intersection(active)
    doc={"start":start,"finish":finish,"parent":pid,"wall_s":time.monotonic()-wall,
         "worker_exit_verified":True,"shards":out,"samples":samples,"radius":R,"seed":SEED,"mp_dps":DPS,
         "python":sys.version,"platform":platform.platform(),"numpy":np.__version__,
         "scipy":scipy_version,"mpmath":mp.__version__}
    Path(a.out).write_text(json.dumps(doc,indent=2)+"\n")
    print(json.dumps({"finish":finish,"wall_s":doc["wall_s"],"workers_exited":True}),flush=True)


if __name__=="__main__":
    main()
