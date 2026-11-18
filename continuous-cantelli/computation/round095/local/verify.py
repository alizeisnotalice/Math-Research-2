"""Independent Mills or adaptive quadrature for the reflected wall convolution."""
import json
import math
from pathlib import Path
import struct
import sys

import mpmath
import scipy
from scipy.integrate import quad
from scipy.special import erfcx, log_ndtr


def f32(x):
    return struct.unpack("f",struct.pack("f",x))[0]


def log_D(z):
    if z>=0:
        return .5*math.log(math.pi/2)+math.log(erfcx(z/math.sqrt(2)))
    return .5*math.log(2*math.pi)+z*z/2+log_ndtr(-z)


def analytic(k,v):
    if v==0:
        return math.log1p(-k*math.exp(log_D(k)))
    lo,hi=log_D(k-v),log_D(k+v)
    return lo+math.log(-math.expm1(hi-lo))-math.log(2*v)


def quad_log(k,v):
    a=k-v
    peak=max(0,-a)
    d=.5*peak*peak

    def integrand(w):
        if w<=0:
            return 0.0
        z=v*w
        if z>24:
            exponent=-.5*(w-peak)**2-max(a,0)*w-math.log(2*v)+math.log1p(-math.exp(-2*z))
        else:
            if abs(z)<.3:
                sinch=sum(z**(2*j)/math.factorial(2*j+1) for j in range(10))
            else:
                sinch=math.sinh(z)/z
            exponent=math.log(w)-w*w/2-k*w-d+math.log(sinch)
        return math.exp(exponent)

    if a>=0:
        scale=1/(1+a)
        value,err=quad(lambda u:scale*integrand(u*scale),0,math.inf,epsabs=1e-14,epsrel=3e-13,limit=200)
    else:
        pos,e1=quad(lambda y:integrand(peak+y),0,math.inf,epsabs=1e-14,epsrel=3e-13,limit=200)
        neg,e2=quad(lambda y:integrand(peak-y),0,peak,epsabs=1e-14,epsrel=3e-13,limit=200)
        value,err=pos+neg,e1+e2
    return math.log(value)+d,err


def precise(k0,v0):
    mp=mpmath.mp
    mp.dps=70
    k,v=map(mp.mpf,(k0,v0))
    def D(z):
        return mp.sqrt(mp.pi/2)*mp.exp(z*z/2)*mp.erfc(z/mp.sqrt(2))
    R=1-k*D(k) if v==0 else (D(k-v)-D(k+v))/(2*v)
    return {"logR":mp.nstr(mp.log(R),55),"k2R":mp.nstr(k*k*R,55)}


def verify(hit):
    k,v=(f32(hit["p"][i]) for i in ("k","v"))
    q,err=quad_log(k,v)
    q0,_=quad_log(k,0)
    D=analytic(k,v)
    gpu=f32(hit["logR"])
    cpu={"logR":q,"logL":q0,"old_factor":math.exp(min(700,q+2*math.log(k))),"estimated_scaled_quad_error":err}
    if k>v:
        upper=q+math.log(k-v)+math.log(k+v)
        err_bound=1/(k-v)**2 + 1/(k*k-v*v) + 1/(k+v)**2
        cpu.update(rational_upper_log_gap=upper,
                   rational_lower_valid=err_bound>=1 or upper>=math.log1p(-err_bound)-2e-12)
    return {"input":{"k":k,"v":v},"gpu_logR":gpu,"gpu_logL":f32(hit["logLinear"]),
            "cpu":cpu,"cpu_mills_logR":D,"gpu_logR_error":abs(gpu-q),"cpu_mills_error":abs(D-q)}


def physical_controls():
    rows=[]
    for name,tau,b,mode in (("fixed",1.,1.,"constant"),("shallow",1.,1.,"moving")):
        values=[]
        for h in (1e-2,1e-4,1e-6,1e-8):
            d=.5 if mode=="constant" else math.sqrt(h)
            sig=math.sqrt(tau*h/(tau+h))
            k,v=d*sig/h,b*sig/tau
            r,_=quad_log(k,v)
            values.append({"h":h,"d":d,"k":k,"v":v,"k2R":math.exp(r+2*math.log(k))})
        rows.append({"case":name,"values":values})
    direct=[]
    for tau,h,b,d in ((1,.03,.8,.25),(.3,.7,1.2,.18),(2,.18,.031,1.25)):
        sigma=math.sqrt(tau*h/(tau+h))
        k,v=d*sigma/h,b*sigma/tau
        y=b+d
        density=lambda variance,z:math.exp(-z*z/(2*variance))/math.sqrt(2*math.pi*variance)
        expected=quad(lambda x:density(h,y-x)*(density(tau,x)-density(tau,2*b-x)),
                      -math.inf,b,epsabs=1e-15,epsrel=2e-12)[0]
        log_r,_=quad_log(k,v)
        normalized=2*v*sigma*density(h,d)*density(tau,b)*math.exp(log_r)
        direct.append({"tau":tau,"h":h,"b":b,"distance":d,
                       "direct_heat":expected,"normalized_heat":normalized,
                       "relative_difference":abs(normalized/expected-1)})
    return {"sequences":rows,"finite_k64_v63":precise(64,63),"direct_physical_checks":direct}


def main():
    raw=json.loads(Path(sys.argv[1]).read_text())
    out={"source":Path(sys.argv[1]).name,"python":sys.version,"scipy":scipy.__version__,
         "mpmath":mpmath.__version__,
         "maxima":{k:verify(v["best"]) for k,v in raw["summaries"].items() if v["best"]},
         "controls":physical_controls()}
    for k in ("naive_bound_failures","small_wall","transition"):
        if k in out["maxima"]:
            row=out["maxima"][k]
            row["mpmath_70digit_working_precision"]=precise(row["input"]["k"],row["input"]["v"])
    Path(sys.argv[2]).write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")


if __name__=="__main__":
    main()
