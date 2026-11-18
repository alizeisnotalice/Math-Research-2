"""Independent displaced-Gaussian quadrature and small-shift checks."""
import json
import math
from pathlib import Path
import struct
import sys
import warnings

import mpmath
import numpy as np
import scipy
from scipy.integrate import IntegrationWarning, quad
from scipy.optimize import brentq


def f32(x):
    return struct.unpack("f",struct.pack("f",x))[0]


def normalize(p):
    return {**p,**{k:f32(p[k]) for k in ("beta","r","g","amp","t")}}


def poly(r,g):
    return np.polynomial.Polynomial([1-3*r,3*g,6*r-1,-g,-r])


def stationary(p,g):
    b=p["beta"]
    P=poly(p["r"],g)
    derivative=P.deriv()-b/(1+b)*np.polynomial.Polynomial([0,1])*P
    return sorted(float(y.real)/math.sqrt(1+b) for y in derivative.roots() if abs(y.imag)<1e-10)


def log_phi2(z):
    if abs(z)<.125:
        return math.log(sum(z**k/math.factorial(k+2) for k in range(14)))
    if z>=12:
        return z+math.log1p(-(1+z)*math.exp(-z))-2*math.log(z)
    return math.log(math.expm1(z)-z)-2*math.log(abs(z))


def coefficients(p,g):
    beta,r,A,t=(p[k] for k in ("beta","r","amp","t"))
    B=1+beta
    xi=t/math.sqrt(B)
    s=t*t/2
    m=math.log(2*A)-1.5*math.log(B)-beta*t*t/(2*B)+math.log(1+g*xi+r*xi*xi)
    W=1/math.sqrt(B)
    roots=stationary(p,g)
    slots=[t+k for k in (-12,-6,-3,-1,0,1,3,6,12)]+[0]
    slots+= [x+k*W/max(1,t) for x in roots+[0] for k in (-10,-4,-1,0,1,4,10)]
    cuts=sorted(set([t-28,t+28]+[x for x in slots if t-28<x<t+28]))
    P=poly(r,g)

    def fun(x):
        y=x/W
        v=float(P(y))
        if v==0:
            return 0.0
        b=math.log(A)+math.log(abs(v))-beta*x*x/2
        z=math.copysign(s*math.exp(b),v)
        return math.exp(-.5*(x-t)**2-.5*math.log(2*math.pi)+2*b+log_phi2(z)-m)

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always",IntegrationWarning)
        ratio,err=quad(fun,cuts[0],cuts[-1],points=cuts[1:-1],epsabs=2e-13,epsrel=3e-12,limit=700)
    if w:
        raise RuntimeError(str(w[0].message))
    diff=math.log(ratio)
    sup=max([0]+[A*math.exp(-beta*x*x/2)*float(P(x/W)) for x in roots])
    return {"log_ratio":diff,"log_scaled_J":diff+m,"log_scaled_mean":m,
            "ratio":ratio,"quad_estimated_absolute_error":err,"b_sup_numeric":sup}


def high_precision(p,g):
    mp=mpmath.mp
    mp.dps=70
    beta,r,A,t,g=map(mp.mpf,(p["beta"],p["r"],p["amp"],p["t"],g))
    B=1+beta
    xi=t/mp.sqrt(B)
    s=t*t/2
    logm=mp.log(2*A)-mp.mpf("1.5")*mp.log(B)-beta*t*t/(2*B)+mp.log(1+g*xi+r*xi*xi)
    cuts=sorted({t-28,t+28,*[t+j for j in (-9,-3,-1,0,1,3,9)],
                  *[k/mp.sqrt(B) for k in (-12,-6,-2,-1,0,1,2,6,12)]})
    cuts=[x for x in cuts if t-28<=x<=t+28]

    def f(x):
        y=mp.sqrt(B)*x
        u=A*mp.exp(-beta*x*x/2)*(1-y*y-g*y*(y*y-3)-r*(y**4-6*y*y+3))
        z=s*u
        if abs(z)<mp.mpf("1e-8"):
            phi=sum(z**k/mp.factorial(k+2) for k in range(12))
        else:
            phi=(mp.expm1(z)-z)/z**2
        return mp.exp(-(x-t)**2/2-logm)*u*u*phi/mp.sqrt(2*mp.pi)

    ratio=mp.quad(f,cuts)
    return {"ratio":mp.nstr(ratio,55),"log_ratio":mp.nstr(mp.log(ratio),55)}


def check(row):
    p=normalize(row["p"])
    pos=coefficients(p,p["g"])
    neg=coefficients(p,-p["g"])
    gpos=f32(row["positive"])
    gneg=f32(row["negative"])
    return {"point":p,"index":row["index"],"gpu_log_ratios":[gpos,gneg],
          "double":[pos,neg],"absolute_log_errors":[abs(gpos-pos["log_ratio"]),abs(gneg-neg["log_ratio"])]}


def controls():
    p={"beta":8.,"r":0.,"g":0.,"amp":.5,"t":.01,"family":0}
    def f(t):
        return coefficients(dict(p,t=t),0)["log_ratio"]
    root=brentq(f,.01,12.)
    return {
        "beta8_Ahalf_balance":root,
        "lower_bracket":high_precision(dict(p,t=root-.001),0),
        "upper_bracket":high_precision(dict(p,t=root+.001),0),
        "ratios_for_global_undercompensation_A1_64":[
            {"t":t,"ratio":coefficients(dict(p,t=t,amp=1/64),0)["ratio"]} for t in (.01,1,3,6,12,24)
        ],
    }


def main():
    gpu=json.loads(Path(sys.argv[1]).read_text())
    out={"input":Path(sys.argv[1]).name,"python":sys.version,"scipy":scipy.__version__,
         "mpmath":mpmath.__version__,
         "maxima":{k:check(v["best"]) for k,v in gpu["summaries"].items() if v["best"]},
         "controlled":controls()}
    if "--precise" in sys.argv:
        for k in ("near_balance","skew_opposite_signs"):
            if k not in out["maxima"]:
                continue
            q=out["maxima"][k]
            p=q["point"]
            q["mpmath_70digit_working_precision"]=[high_precision(p,g) for g in (p["g"],-p["g"])]
    Path(sys.argv[2]).write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")


if __name__=="__main__":
    main()
