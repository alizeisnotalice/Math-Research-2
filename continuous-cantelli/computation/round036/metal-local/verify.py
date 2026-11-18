"""Independent adaptive integration of the GPU maxima and exact finite witnesses."""
from dataclasses import dataclass
import json
import math
from pathlib import Path
import struct
import sys
import warnings

import mpmath
import scipy
from scipy.integrate import IntegrationWarning, quad


@dataclass
class Point:
    m: float
    a: float
    r: float
    base_edge: float
    shape: float
    power: float
    epsilon: float
    delta: float
    left: float
    right: float

    def q(self, side, baseline):
        return 0 if baseline else self.left if side < 0 else self.right


def variance(p, u, q):
    if u <= 0:
        return p.epsilon
    t = u/p.r
    z = p.power*max(0.0, math.log(1/t))
    factor = 1.0 if q == 0 or z > math.log1p(750/p.shape) else 1+q*math.exp(-p.shape*math.expm1(z))
    return p.epsilon+p.base_edge*t**4*factor


def integrate(p, baseline=False):
    target = p.m+p.epsilon
    y = p.a+p.delta
    offset = 0.5*y*y/target
    pref = math.sqrt(target/(2*math.pi))
    result, errors = [], []
    for side in (-1,1):
        q = p.q(side, baseline)
        center = side*p.delta
        lo, hi = -center, p.r-center
        cuts = {lo,hi}
        for scale in (math.sqrt(variance(p,min(p.r,max(0,center)),q)), p.r*p.r/math.sqrt(p.base_edge),
                      p.r*(p.epsilon/p.base_edge)**.25):
            for factor in (-32,-16,-8,-4,-2,-1,0,1,2,4,8,16,32):
                w = factor*scale
                if lo < w < hi:
                    cuts.add(w)
                w = factor*scale-center
                if lo < w < hi:
                    cuts.add(w)
        for drop in (1/8,1/4,1/2,1,2,3,4,6,10,20,30,50):
            w = p.r*math.exp(-math.log1p(drop/p.shape)/p.power)-center
            if lo < w < hi:
                cuts.add(w)

        def f(w):
            u = max(0.0,min(p.r,center+w))
            v = variance(p,u,q)
            x = p.a+side*u
            z = w/math.sqrt(v)
            return pref*math.exp(offset-(x*x+z*z)/2)/math.sqrt(v)

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", IntegrationWarning)
            value, estimated = quad(f,lo,hi,points=sorted(cuts)[1:-1],
                                    epsabs=1e-11,epsrel=2e-12,limit=500)
        if caught:
            raise RuntimeError(f"Adaptive quad warning {caught[0].message} at {p}")
        result.append(value)
        errors.append(estimated)
    return {"ratio":sum(result),"sides":result,"estimated_errors":errors}


def high_precision(p, baseline=False):
    mp = mpmath.mp
    mp.dps = 65
    m,a,r,base_edge,shape,power,eps,delta,left,right = map(mp.mpf,vars(p).values())
    target = m+eps
    offset = (a+delta)**2/(2*target)
    pref = mp.sqrt(target/(2*mp.pi))
    result = []
    for side in (-1,1):
        q = mp.mpf(0) if baseline else (left if side < 0 else right)
        center = side*delta
        lo,hi = -center,r-center
        cuts = {lo,hi,0}

        def variance(u):
            if u <= 0:
                return eps
            t = u/r
            return eps+base_edge*t**4*(1+q*mp.exp(shape*(1-t**(-power))))

        for scale in (mp.sqrt(variance(max(0,min(r,center)))), r*r/mp.sqrt(base_edge), r*(eps/base_edge)**(mp.mpf(1)/4)):
            for factor in (-32,-16,-8,-4,-2,-1,0,1,2,4,8,16,32):
                cuts.add(factor*scale)
                cuts.add(factor*scale-center)
        for drop in (1/4,1/2,1,2,4,8,16,32):
            cuts.add(r*(1+drop/shape)**(-1/power)-center)

        def f(w):
            u = max(0,min(r,center+w))
            x = a+side*u
            v = variance(u)
            return pref*mp.exp(offset-(x*x+w*w/v)/2)/mp.sqrt(v)

        result.append(mp.quad(f,sorted(w for w in cuts if lo<=w<=hi)))
    return mp.nstr(sum(result),59)


def f32(value):
    return struct.unpack("f", struct.pack("f", value))[0]


def verify_hit(hit):
    raw = hit["point"]
    point = Point(*(f32(raw[k]) for k in ("m","a","r","baseEdge","shape","power","epsilon","delta","left","right")))
    out = {}
    threshold = point.m*math.log(point.m)/(point.m-1)
    for k, base in (("ratio",False),("baseline",True)):
        result = integrate(point,base)
        out[k] = {
            "gpu":f32(hit[k]),"cpu":result,"relative_error":abs(f32(hit[k])/result["ratio"]-1),
        }
    out["true_margin"] = point.a*point.a-threshold
    out["admissible"] = point.base_edge*max(1,1+point.left,1+point.right) <= point.m
    return out


def main():
    gpu = json.loads(Path(sys.argv[1]).read_text())
    result = {
        "input":Path(sys.argv[1]).name,"python":sys.version,"scipy":scipy.__version__,"mpmath":mpmath.__version__,
        "maxima":{name:verify_hit(row["best"]) for name,row in gpu["summaries"].items()},
    }
    Path(sys.argv[2]).write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")


if __name__ == "__main__":
    main()
