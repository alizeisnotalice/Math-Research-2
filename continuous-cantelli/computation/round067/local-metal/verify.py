"""Independent adaptive quadrature for the periodic heat scan."""
import cmath
import json
import math
from pathlib import Path
import struct
import sys
import warnings

import mpmath
import scipy
from scipy.integrate import IntegrationWarning, quad


F32_KEYS = ("c1","s1","c2","s2","c3","s3","b1","t1","b2","t2","b3","t3","d","eta","p")


def f32(x):
    return struct.unpack("f", struct.pack("f", x))[0]


def read_point(p):
    return {**p, **{k:f32(p[k]) for k in F32_KEYS}}


def trig(x, coefficients):
    return sum(c*math.cos(i*x)+s*math.sin(i*x)
               for i,(c,s) in enumerate(zip(coefficients[::2], coefficients[1::2]),1))


def amplitudes(p, mp=None):
    c = [p[k] for k in ("c1","s1","c2","s2","c3","s3")]
    b = [p[k] for k in ("b1","t1","b2","t2","b3","t3")]
    return (c,b) if mp is None else (list(map(mp.mpf,c)),list(map(mp.mpf,b)))


def coefficients(p, force_exact_zero=False):
    g,n,eta,d,power = (p[key] for key in ("g","n","eta","d","p"))
    scaled = n%g != 0 and eta>0
    if force_exact_zero and n%g and eta==0:
        return 0j,0j
    c,b = amplitudes(p)
    width = (2/(d*n*n))**(1/power)
    end = min(math.pi/g, (math.pi/2)*width*(61/.35)**(1/power))
    cuts = {0.0,end}
    cuts.update(val for exponent in (-8,-4,-2,0,1,2,3,4,5,6,7)
                if 0<(val:=width*2**(exponent/power))<end)
    pieces = math.ceil(n*end/10)
    cuts.update(end*i/pieces for i in range(1,pieces))
    cuts=sorted(cuts)
    absolute = min(2e-13,width*1e-9)
    total = [0.0]*4

    for j in range(g):
        origin = 2*math.pi*j/g
        phase = 2*math.pi*(n*j%g)/g
        for sign in (-1,1):
            def integrand(x,k):
                rate = d*n*n/2*(2*math.sin(g*x/2)/g)**power
                u = sign*x
                q = 1+trig(g*u,c)
                h = trig(origin+u,b)
                base = math.exp(-rate*q)
                if scaled:
                    shift = -rate*eta*h
                    if abs(shift)<1:
                        target = base*math.expm1(shift)/eta
                    else:
                        target = (math.exp(-rate*(q+eta*h))-base)/eta
                    linear = -rate*h*base
                else:
                    target,linear = math.exp(-rate*(q+eta*h)),0
                arg = phase+n*u
                return (math.cos(arg) if k%2==0 else math.sin(arg))*(target if k<2 else linear)

            for k in range(4 if scaled else 2):
                with warnings.catch_warnings(record=True) as record:
                    warnings.simplefilter("always", IntegrationWarning)
                    value,_=quad(integrand,0,end,args=(k,),points=cuts[1:-1],
                                 epsabs=absolute,epsrel=4e-12,limit=2000)
                if record:
                    raise RuntimeError(f"adaptive quadrature warning {record[0].message}")
                total[k]+=value/(2*math.pi)

    return complex(*total[:2]),complex(*total[2:])


def precise(p):
    mp=mpmath.mp
    mp.dps=65
    n,g=p["n"],p["g"]
    eta,d,power=(mp.mpf(p[key]) for key in ("eta","d","p"))
    c,b=amplitudes(p,mp)
    s=(2/(d*n*n))**(1/power)
    end=min(mp.pi/g,mp.pi*s*(mp.mpf(220)/mp.mpf(".35"))**(1/power)/2)
    limits=sorted({mp.mpf(0),end,*(z for k in (-8,-4,-2,0,1,2,3,4,5,6,7)
                                  if 0<(z:=s*mp.power(2,mp.mpf(k)/power))<end)})
    total=mp.mpc(0)
    for j in range(g):
        origin=2*mp.pi*j/g
        phase=2*mp.pi*(n*j%g)/g
        for side in (-1,1):
            def fun(x):
                rate=d*n*n/2*(2*mp.sin(g*x/2)/g)**power
                u=side*x
                q=1+sum(c[2*k-2]*mp.cos(k*g*u)+c[2*k-1]*mp.sin(k*g*u) for k in (1,2,3))
                h=sum(b[2*k-2]*mp.cos(k*(origin+u))+b[2*k-1]*mp.sin(k*(origin+u)) for k in (1,2,3))
                weight=mp.exp(-rate*q)*mp.expm1(-rate*eta*h)/eta if n%g and eta else mp.exp(-rate*(q+eta*h))
                return weight*mp.exp(1j*(phase+n*u))

            total+=mp.quad(fun,limits)/(2*mp.pi)
    return [mp.nstr(total.real,54),mp.nstr(total.imag,54)]


def vec(z):
    return [z.real,z.imag]


def verify_hit(hit):
    point=read_point(hit["p"])
    exact,linear=coefficients(point)
    gpu=complex(f32(hit["real"]),f32(hit["imag"]))
    gpu_l=complex(f32(hit["linearReal"]),f32(hit["linearImag"]))
    return {
        "point":point,"gpu":vec(gpu),"gpu_linear":vec(gpu_l),"cpu":vec(exact),"cpu_linear":vec(linear),
        "absolute_error":abs(gpu-exact),"linear_absolute_error":abs(gpu_l-linear),
    }


def asymptotic_constant(p,n):
    g,d,power,eta=(p[key] for key in ("g","d","p","eta"))
    c,b=amplitudes(p)
    total=0j
    for j in range(g):
        origin=2*math.pi*j/g
        q=1+trig(0,c)+eta*trig(origin,b)
        if power<2:
            weight=math.gamma(1+1/power)*(2/(d*q))**(1/power)/math.pi
        elif power==2:
            weight=math.exp(-1/(2*d*q))/math.sqrt(2*math.pi*d*q)
        else:
            weight=0
        total+=cmath.exp(2j*math.pi*(n*j%g)/g)*weight
    if n%g and eta:
        total/=eta
    return total


def controlled():
    base=dict(zip(F32_KEYS, (.12,.05,-.08,.04,.06,-.03,.5,-.2,.1,.1,-.06,.04,1.0,.125,1.0)))
    results=[]
    for g,power in ((1,.75),(1,2),(1,4),(3,1),(3,2)):
        p=dict(base,g=g,p=power)
        records=[]
        for n in (25,97,385):
            p["n"]=n
            f,_=coefficients(p)
            actual=n**(2/power)*f
            records.append({
                "n":n,"scaled":vec(actual),"limit":vec(asymptotic_constant(p,n)),
                "absolute_difference":abs(actual-asymptotic_constant(p,n)),
            })
        results.append({"g":g,"p":power,"profile":dict(p),"values":records})
    return results


def main():
    ledgers=json.loads(Path(sys.argv[1]).read_text())
    data={"source":Path(sys.argv[1]).name,"python":sys.version,
          "scipy":scipy.__version__,"mpmath":mpmath.__version__,
          "maxima":{k:verify_hit(v["best"]) for k,v in ledgers["summaries"].items() if v["best"]},
          "controlled":controlled()}
    if "--precise" in sys.argv:
        for name in ("isolated_subquadratic_asymptotic","weakly_broken_linear_check"):
            target=data["maxima"][name]
            target["mpmath_65digits"]=precise(target["point"])
    Path(sys.argv[2]).write_text(json.dumps(data,indent=2,sort_keys=True)+"\n")


if __name__=="__main__":
    main()
