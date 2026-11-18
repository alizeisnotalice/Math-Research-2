"""Independent integration of Metal representatives and exactly tuned profiles."""
import json
import math
from pathlib import Path
import struct
import sys

import mpmath
import scipy
from scipy.integrate import quad


def f32(x):
    return struct.unpack("f",struct.pack("f",x))[0]


def knots(p):
    r, edge, eps=(p[key] for key in ("r","edge","eps"))
    points={0.0,r}
    for scale in (math.sqrt(eps),r*r/math.sqrt(edge),r*(eps/edge)**.25):
        points.update(t for factor in (1/32,1/8,1/2,1,2,4,8,16,32)
                      if 0<(t:=scale*factor)<r)
    return sorted(points)


def integrate(p):
    r,m,edge,b,eps,h,k=(p[key] for key in ("r","m","edge","b","eps","h","k"))
    a=k*h/r
    target=m+eps
    normal=2*math.sqrt(target/(2*math.pi))*math.exp(-(a*a/2)*(1-1/target))

    def integrand(u,index):
        t=u/r
        c=edge*t**4
        v=eps+c*(1+b*t*t)
        z=u*u/v
        if z>800:
            return 0.0
        slope=c*t/v
        base=math.exp(-(u*u+z)/2)/math.sqrt(v)
        lead=(z-1)*slope/2
        if index==0:
            return base*(z*z-6*z+3)*slope*slope/8
        if index==1:
            return base*t*lead
        if index==2:
            x=h*slope
            xx=x*x
            even=-math.log1p(-xx)/4-z*xx/(2*(1-xx))
            odd=(z*x/(1-xx)-math.atanh(x))/2
            au=k*h*t
            result=(math.cosh(au)*(math.expm1(even)/h**2
                        +2*math.exp(even)*(math.sinh(odd/2)/h)**2)
                    -(math.sinh(au)/h)*math.exp(even)*math.sinh(odd)/h)
            return base*result
        d=k*t
        e2=(.25-z/2)*slope*slope
        e4=(.125-z/2)*slope**4
        q1=d-lead
        q3=-(z-1/3)*slope**3/2
        return base*(e4+e2*e2/2+e2*q1*q1/2+q1*q3+(q1**4-d**4)/24)

    cuts=knots(p)[1:-1]
    vals=[quad(integrand,0,r,args=(j,),points=cuts,epsabs=2e-12,epsrel=2e-12,limit=300) for j in range(4)]
    A,B,X,J4=(v[0] for v in vals)
    return {"A":A,"B":B,"root":A/B if B else None,
            "exact":normal*X,"leading":normal*(A-k*B),"J4":J4,
            "estimated_errors":[v[1] for v in vals]}


def verify_hit(hit):
    p={k:(v if k=="tuned" else f32(v)) for k,v in hit["point"].items()}
    q=integrate(p)
    return {"point":p,"cpu":q,"gpu_exact":hit["baseline"],
            "gpu_gap":hit["ratio"],
            "exact_absolute_error":abs(q["exact"]-f32(hit["baseline"])),
            "gap_absolute_error":abs(abs(q["exact"]-q["leading"])-f32(hit["ratio"]))}


def mp_check(p, hs):
    mp=mpmath.mp
    mp.dps=65
    r,m,edge,b,eps=(mp.mpf(p[k]) for k in ("r","m","edge","b","eps"))
    cuts=list(map(mp.mpf,knots(p)))

    def params(u):
        t=u/r
        c=edge*t**4
        v=eps+c*(1+b*t*t)
        z=u*u/v
        slope=c*t/v
        base=mp.exp(-(u*u+z)/2)/mp.sqrt(v)
        return t,c,v,z,slope,base

    def coeff(u,index):
        t,c,v,z,slope,base=params(u)
        return base*((z*z-6*z+3)*slope*slope/8 if index==0 else t*(z-1)*slope/2)

    A,B=[mp.quad(lambda u:coeff(u,j),cuts) for j in range(2)]
    k=A/B
    T=m+eps
    norm0=2*mp.sqrt(T/(2*mp.pi))

    def pointwise(u,h):
        t,c,v,z,slope,base=params(u)
        au=k*h*t
        summed=mp.mpf(0)
        for sign in (-1,1):
            new_v=v+sign*h*c*t
            new=mp.exp(-(u*u+u*u/new_v)/2)/mp.sqrt(new_v)
            summed += mp.exp(-sign*au)*(new-base)
        return summed

    def j4(u):
        t,c,v,z,slope,base=params(u)
        d=k*t
        e2=(mp.mpf(1)/4-z/2)*slope*slope
        e4=(mp.mpf(1)/8-z/2)*slope**4
        q1=d-(z-1)*slope/2
        q3=-(z-mp.mpf(1)/3)*slope**3/2
        return base*(e4+e2*e2/2+e2*q1*q1/2+q1*q3+(q1**4-d**4)/24)

    quarter=norm0*mp.quad(j4,cuts)
    observed={}
    for h0 in hs:
        h=mp.mpf(h0)
        a=k*h/r
        norm=mp.sqrt(T/(2*mp.pi))*mp.exp(-a*a*(1-1/T)/2)
        observed[str(h0)]=mp.nstr(norm*mp.quad(lambda u:pointwise(u,h),cuts)/h**4,40)
    return {"k_star":mp.nstr(k,40),"fourth_coefficient":mp.nstr(quarter,40),"ratio_over_h4":observed}


def main():
    gpu=json.loads(Path(sys.argv[1]).read_text())
    # Choose controlled positive and negative fourth terms from a small fixed grid.
    choices=[]
    for b in (-.25,0,.5):
        for factor in (.125,.25,.5,1,2,4,8,16):
            p={"m":2.0,"r":.125,"edge":factor*.125**2,"b":b,"eps":2**-20,"h":.01,"k":0.0}
            k=integrate(p)["root"]
            if k is None or abs(k)>12: continue
            p["k"]=k
            choices.append((integrate(p)["J4"],p))
    best=(min(choices,key=lambda x:x[0]),max(choices,key=lambda x:x[0]))
    report={"python":sys.version,"scipy":scipy.__version__,"mpmath":mpmath.__version__,
            "gpu_source":Path(sys.argv[1]).name,
            "maxima":{key:verify_hit(row["best"]) for key,row in gpu["summaries"].items()},
            "controlled":[{"input":p,"precise":mp_check(p,("0.01","0.001","0.0001"))} for _,p in best]}
    Path(sys.argv[2]).write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")


if __name__=="__main__":
    main()
