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

R=.1875
DPS=87
SEED=581119
ORDERS=(4,6,8,10,12,14)
PANELS=(0,1,3,6,10,15,20)
CPU_RE=re.compile(r"CPU usage: (.*?)% user, (.*?)% sys, (.*?)% idle")


def utc():
    return datetime.now(timezone.utc).isoformat(timespec="milliseconds")


def grid(n):
    rows=qmc.Sobol(6,scramble=True,seed=SEED).random_base2((n-1).bit_length())[:n]
    for j,(center,c_s,b_s,b_other,d_s,d_other) in enumerate(rows):
        order=ORDERS[j%len(ORDERS)]
        match=(j//len(ORDERS))%2==0
        a=float(7*center-3.5)
        c=float(.3+2.35*c_s)
        b=float(1.4*b_s-.7)
        b_minus=-b if match else float(1.4*b_other-.7)
        yield j,a,(order,order),([c,b,float(.8*d_s-.35)],[c,b_minus,float(.8*d_other-.35)]),match


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
    a,c,b,d=(mp.mpf(q) for q in (a,c,b,d))
    C0=c*n*(n-1)/2
    C1=c*n*(n+1)*(b-sig*a)/2
    C2=c*(n+1)*(n+2)*(d-sig*a*b+(a*a-1)/2)/2
    if n==4:
        C2+=c*c*mp.factorial(2*n)/(8*mp.factorial(2*n-4))
    return C0,C1,C2


def run(rows,smoke):
    mp.mp.dps=DPS
    start=utc(); cpu=time.process_time()
    result=[]
    for idx,a,orders,shape,match in rows:
        n=orders[0]
        exps=(5,7,9,12,15,18) if n<8 else (4,6,8,10,12,14)
        if smoke:
            exps=exps[0],exps[-1]
        info=[]
        C=[leading(n,a,*data,sig) for data,sig in zip(shape,(1,-1))]
        order=n if match else n-1
        fac=mp.factorial(order)
        target=fac*(C[0][2]-C[1][2] if match else C[0][1]+C[1][1])
        for ex in exps:
            z=R*2**(-ex)
            Z=mp.mpf(z)
            halves=[]
            for data,sig in zip(shape,(1,-1)):
                h=high(n,a,*data,z,sig)
                f=fast(n,a,*data,z,sig)
                err=abs(float(h)-f)
                if err>3e-10*max(1,abs(f)):
                    raise ArithmeticError((idx,n,ex,str(h),f,err))
                halves.append((h,err))
            left,right=halves[0][0],halves[1][0]
            if match:
                jj=fac*((left-C[0][0]-Z*C[0][1])-(right-C[1][0]-Z*C[1][1]))/Z**2
            else:
                jj=fac*(left-C[0][0]+right-C[1][0])/Z
            info.append({"exp":ex,"z":z,"observed_jump":mp.nstr(jj,39),
                         "jump_error":mp.nstr(jj-target,30),
                         "sides":[{"sig":sig,"order":n,"normalized":mp.nstr(h,41),
                                   "expect0":mp.nstr(C[j][0],31),"expect1":mp.nstr(C[j][1],31),
                                   "expect2":mp.nstr(C[j][2],31),
                                   "quadratic_error":mp.nstr((h-C[j][0]-Z*C[j][1])/Z**2-C[j][2],26),
                                   "cross_error":err}
                                  for j,((h,err),sig) in enumerate(zip(halves,(1,-1)))]})
        result.append({"index":idx,"center":a,"orders":orders,"coefficients":shape,
                       "linear_terms_matched":match,"jump_derivative_order":order,
                       "predicted_jump_divided_by_source":mp.nstr(target,33),"rows":info})
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
