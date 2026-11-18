#!/usr/bin/env python3
"""Adaptive whole-line CPU audit of selected noncompact profile tests."""
import json
import re
import statistics
import struct
import sys
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import quad

NORM=1/np.sqrt(2*np.pi)


def f32(x):
    return struct.unpack("=f",struct.pack("=f",x))[0]


def tail(x,r,a,p):
    inv=r/x
    return a*inv**p*(1-inv)


def data_for(profile):
    b=np.array(profile["offset"],dtype=np.float32).astype(float)
    dx=profile["dx"]
    x=np.linspace(-profile["R"],profile["R"],len(b))
    p,w=np.polynomial.legendre.leggauss(32)
    p,w=(p+1)/2,w/2
    return ravel(x[:-1,None]+dx*p),ravel(b[:-1,None]+np.diff(b)[:,None]*p),np.tile(dx*NORM*w,len(b)-1)


def ravel(a):
    return a.ravel()


def integrate(function,l,r):
    return quad(function,l,r,epsabs=4e-14,epsrel=1e-12,limit=150)[0]


def moderate(core,z,R,ap,am,pp,pm):
    x,b,w=core
    aa,bb=z.real,z.imag
    def integrand(x,b):
        dz=x-aa
        return np.exp(-dz*dz/2+1j*bb*dz)*np.expm1(-z*z*b/2)
    compact=np.dot(w,integrand(x,b))
    extra=0j
    for sg,A,P in ((1,ap,pp),(-1,am,pm)):
        def fn(r):
            return NORM*integrand(sg*r,tail(r,R,A,P))
        extra+=integrate(lambda r:fn(r).real,R,np.inf)
        extra+=1j*integrate(lambda r:fn(r).imag,R,np.inf)
    return compact+extra,extra


def far(core,t,R,ap,am,pp,pm):
    x,b,w=core
    q=t*t/2
    compact=np.dot(w,np.exp(-(x-t)**2/2-q*b))
    whole=compact
    for sg,A,P in ((1,ap,pp),(-1,am,pm)):
        lo=R-sg*t
        bounds=[lo]+[z for z in (-25,-10,-4,0,4,10,25) if lo<z]+[np.inf]
        vals=[integrate(lambda z:NORM*np.exp(-z*z/2-q*tail(z+sg*t,R,A,P)),l,r)
              for l,r in zip(bounds,bounds[1:])]
        whole+=sum(vals)
    return float(np.log(whole))


def main():
    profiles=json.loads(Path("candidates.json").read_text())
    data=json.loads(Path(sys.argv[1]).read_text())
    records=[]
    for profile,r in zip(profiles,data["byProfile"],strict=True):
        core=data_for(profile)
        modes={}
        for i,(mode,ledger) in enumerate(zip(data["modes"],r,strict=True)):
            h=ledger["best"]
            p,q,a,b,ap,am,pp,pm=map(f32,(h["p"],h["q"],h["a"],h["b"],h["ap"],h["am"],h["pp"],h["pm"]))
            R=profile["R"]
            if i<4:
                measured,extra=moderate(core,complex(p,q),R,ap,am,pp,pm)
                if i<2: assert p not in profile["training"]
                check=dict(cpu=[measured.real,measured.imag],gpu=[a,b],
                           abs_error=abs(measured-(a+1j*b)),tail_only=[extra.real,extra.imag],
                           compact_only=[(measured-extra).real,(measured-extra).imag])
            else:
                measured=far(core,p,R,ap,am,pp,pm)
                power=pp if p>0 else pm
                pref=-.5*(ap if p>0 else am)*R**power
                check=dict(cpu=measured,gpu=a,abs_error=abs(measured-a),
                           tilted_center_log=b,scaled_cpu=measured/abs(p)**(2-power),
                           asymptotic_scaled_log=pref,sign_matches_expected=measured*pref>0)
            modes[mode]=dict(p=p,q=q,ap=ap,am=am,pp=pp,pm=pm,tested=ledger["tested"],
                             negative=ledger["negative"],nonfinite=ledger["nonfinite"],**check)
        records.append(dict(case=profile["case"],modes=modes))
    stats={}
    if len(sys.argv)>3:
        v=[int(m.group(1)) for l in Path(sys.argv[3]).read_text().splitlines()
           if l.startswith("during ") and (m:=re.search(r'"Device Utilization %"=(\d+)',l))]
        stats=dict(count=len(v),median=statistics.median(v),min=min(v),max=max(v))
    errors={mode:max(r["modes"][mode]["abs_error"] for r in records) for mode in data["modes"]}
    out=dict(numpy=np.__version__,scipy=scipy.__version__,records=records,max_error_by_mode=errors,
             ioreg=stats,count=sum(m["tested"] for r in records for m in r["modes"].values()))
    Path(sys.argv[2]).write_text(json.dumps(out,indent=2)+"\n")
    print("errors",errors)
    for row in records:
        print(row["case"],{k:(v["p"],v["q"],v["cpu"]) for k,v in row["modes"].items()})
    print("IOReg",stats)


if __name__=="__main__":
    main()
