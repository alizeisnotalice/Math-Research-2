"""Adaptive double and multiprecision audits of the Cantor-gap decomposition."""
import json
import math
from pathlib import Path
import struct
import sys
import warnings

import mpmath
import scipy
from scipy.integrate import IntegrationWarning, quad
from scipy.special import gamma, gammainc


def f32(x):
    return struct.unpack("f",struct.pack("f",x))[0]


def normalize(p):
    return {**p,"p":f32(p["p"]),"d":f32(p["d"]),
            "lambda":1/3 if p["ternary"] else f32(p["lambda"])}


def level_integral(power, omega, phase, s):
    error=s*s/(2*(2*power+1))+s*omega*omega/((power+1)*(power+2)*(power+3))
    if error<1e-21:
        return -s/(power+1)
    sub=s<=0.25
    width=s**(-1/power)
    end=1.0 if sub else min(1.0,width*65**(1/power))
    cuts={0.0,end}
    cuts.update(t for k in (-4,-2,0,1,2,3,4,5)
                if 0<(t:=width*2**(k/power))<end)
    pieces=math.ceil(omega*end/12)
    cuts.update(end*i/pieces for i in range(1,pieces))
    def f(t):
        phase_factor=math.cos(phase-omega*t)
        exp=-s*t**power
        return phase_factor*(math.expm1(exp) if sub else math.exp(exp))

    # The first cell uses sqrt(t) as coordinate to resolve half-power cusps.
    ts=sorted(cuts)
    t0=ts[1]
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", IntegrationWarning)
        abs_tol=min(3e-14,s*2e-10) if sub else 3e-14
        val=quad(lambda z:f(t0*z*z)*2*t0*z,0,1,epsabs=abs_tol,epsrel=3e-12)[0]
        if len(ts)>2:
            val+=quad(f,t0,end,points=ts[2:-1],epsabs=abs_tol,epsrel=3e-12,limit=1000)[0]
    if caught:
        raise RuntimeError(str(caught[0].message))
    return val if sub else val-math.sin(phase)/omega


def coefficients(p):
    lam,power,d,n=(p[k] for k in ("lambda","p","d","n"))
    ternary=p["ternary"]!=0
    product=1.0
    coarse,fine,mass_delta=0.0,0.0,0.0
    ratio=2*lam**(power+1)
    for j in range(1,161):
        scale=lam**(j-1)
        length=(1-2*lam)*scale
        den=3**j if ternary else 0
        omega=math.pi*n*length
        phase=math.pi*(n%(2*den))/den if ternary else omega
        s=2*math.pi**2*n*n*d*(length/2)**power
        volume=(1-2*lam)*(2*lam)**(j-1)
        value=level_integral(power,omega,phase,s)
        signed=volume*product*(-1 if n%2 else 1)*value
        is_coarse=(n>=den if ternary else n*length>=1) or s>=1
        if is_coarse:
            coarse+=signed
        else:
            fine+=signed
        if s<1e-4:
            mass_delta+=volume*(-s/(power+1)+s*s/(2*(2*power+1))-s**3/(6*(3*power+1)))
        else:
            mass_delta+=volume*(gamma(1+1/power)*gammainc(1/power,s)*s**(-1/power)-1)
        tail=volume*s*ratio/((1-ratio)*(power+1))
        if tail<1e-19:
            break
        phase_product=2*math.pi*(n%den)/den if ternary else math.pi*n*(1-lam)*scale
        product*=math.cos(phase_product)
    return {"full":coarse+fine,"coarse":coarse,"fine":fine,
            "mass":1+mass_delta,"levels":j,"tail":tail}


def precise(p):
    mp=mpmath.mp
    mp.dps=65
    power,d=map(mp.mpf,(p["p"],p["d"]))
    n=p["n"]
    ternary=p["ternary"]!=0
    lam=mp.mpf(1)/3 if ternary else mp.mpf(p["lambda"])
    product=mp.mpf(1)
    coarse=fine=mp.mpf(0)
    ratio=2*lam**(power+1)
    # Use independent direct gap-center integrals until the tail reaches 1e-25.
    for j in range(1,161):
        scale=lam**(j-1)
        length=(1-2*lam)*scale
        den=3**j if ternary else 0
        omega=mp.pi*n*length
        phase=mp.pi*(n%(2*den))/den if ternary else omega
        s=2*mp.pi**2*n*n*d*(length/2)**power
        volume=(1-2*lam)*(2*lam)**(j-1)
        tail=volume*s*ratio/((1-ratio)*(power+1))
        sub=s<=1
        base=s**(-1/power)
        cap=mp.mpf(1) if sub else min(mp.mpf(1),base*mp.mpf(85)**(1/power))
        splits=max(1,int(mp.ceil(omega*cap/20)))
        cuts={mp.mpf(0),cap,*(cap*k/splits for k in range(1,splits))}
        cuts.update(t for el in (-4,-2,0,1,2,3,4,5,6)
                    if 0<(t:=base*mp.power(2,mp.mpf(el)/power))<cap)
        v=mp.quad(lambda t:mp.cos(phase-omega*t)*
                  (mp.expm1(-s*t**power) if sub else mp.exp(-s*t**power)),sorted(cuts))
        value=v if sub else v-mp.sin(phase)/omega
        signed=volume*product*(-1 if n%2 else 1)*value
        is_coarse=(n>=den if ternary else n*length>=1) or s>=1
        if is_coarse:
            coarse+=signed
        else:
            fine+=signed
        if tail<mp.mpf("1e-25"):
            break
        product*=mp.cos(2*mp.pi*(n%den)/den if ternary else mp.pi*n*(1-lam)*scale)
    return {key:mp.nstr(value,42) for key,value in [("full",coarse+fine),("coarse",coarse),
             ("fine",fine),("rigorous_series_tail_bound",tail)]}


def verify_hit(hit):
    p=normalize(hit["p"])
    actual=coefficients(p)
    gpu={k:f32(hit[k]) for k in ("full","coarse","fine")}
    return {
        "point":p,"gpu":gpu,"cpu":actual, "gpu_depth":hit["levels"],
        "gpu_tail_estimate":f32(hit["tail"]),
        "errors":{k:abs(gpu[k]-actual[k]) for k in gpu},
    }


def cantor_measure_fourier(n):
    return (-1 if n%2 else 1)*math.prod(math.cos(2*math.pi*(n%(3**j))/(3**j)) for j in range(1,85))


def controlled():
    cases=[]
    for power in (.75,1.5,2.0,3.0):
        for shift in (-1,0,1):
            records=[]
            for k in (2,4,7,9):
                p={"lambda":1/3,"p":power,"d":.3,"n":3**k+shift,"ternary":1,"frequency":0}
                result=coefficients(p)
                records.append({"n":p["n"],"k":k,"full":result["full"],"mass":result["mass"],
                                "normalized":result["full"]/result["mass"]})
            cases.append({"p":power,"d":.3,"shift":shift,"records":records,
                          "limit":cantor_measure_fourier(1)*cantor_measure_fourier(shift) if power<2 else None})
    return cases


def main():
    gpu=json.loads(Path(sys.argv[1]).read_text())
    data={"source":Path(sys.argv[1]).name,"python":sys.version,"scipy":scipy.__version__,
          "mpmath":mpmath.__version__,"maxima":{name:verify_hit(row["best"]) for name,row in gpu["summaries"].items() if row["best"]},
          "controlled":controlled()}
    if "--precise" in sys.argv:
        for name in ("ternary_resonant","coarse_fine_cancellation"):
            x=data["maxima"][name]
            x["mpmath_65digit_working_precision"]=precise(x["point"])
    Path(sys.argv[2]).write_text(json.dumps(data,indent=2,sort_keys=True)+"\n")


if __name__=="__main__":
    main()
