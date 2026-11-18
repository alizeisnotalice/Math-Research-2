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
from scipy.special import ndtr


@dataclass
class Point:
    m: float
    a: float
    r: float
    edge: float
    shape: float
    power: float
    epsilon: float
    delta: float
    family: int

    def shape_on(self, side):
        if side < 0:
            return self.shape
        return (self.shape, 0, 2*self.shape)[self.family]


def variance(p, u, shape):
    if shape == 0 or u <= 0:
        return p.epsilon
    log_power = p.power*max(0.0, math.log(p.r/u))
    if log_power > math.log1p(750/shape):
        return p.epsilon
    return p.epsilon+p.edge*math.exp(-shape*math.expm1(log_power))


def integrate(p):
    target = p.m+p.epsilon
    y = p.a+p.delta
    offset = 0.5*y*y/target
    pref = math.sqrt(target/(2*math.pi))
    result, errors = [], []
    for side in (-1, 1):
        shape = p.shape_on(side)
        if shape == 0:
            # Product of the two Gaussians, integrated over the right plateau.
            total = 1+p.epsilon
            posterior_mean_relative = (p.delta-p.a*p.epsilon)/total
            sigma = math.sqrt(p.epsilon/total)
            zlo = -posterior_mean_relative/sigma
            zhi = (p.r-posterior_mean_relative)/sigma
            if zlo > 0:
                probability = ndtr(-zlo)-ndtr(-zhi)
            else:
                probability = ndtr(zhi)-ndtr(zlo)
            result.append(math.sqrt(target/total)*math.exp(offset-y*y/(2*total))*probability)
            errors.append(0.0)
            continue
        center = side*p.delta
        lo, hi = -center, p.r-center
        cuts = {lo, hi}
        radius_scale = math.sqrt(variance(p, min(p.r, max(0,center)), shape))
        for factor in (-32,-16,-8,-4,-2,-1,0,1,2,4,8,16,32):
            w = factor*radius_scale
            if lo < w < hi:
                cuts.add(w)
        max_drop = max(4, math.ceil(math.log(p.edge/p.epsilon))+4)
        for k in range(1, max_drop+1):
            u = p.r*math.exp(-math.log1p(k/shape)/p.power)
            w = u-center
            if lo < w < hi:
                cuts.add(w)

        def f(w):
            u = max(0.0, min(p.r, center+w))
            x = p.a+side*u
            v = variance(p, u, shape)
            z = w/math.sqrt(v)
            return pref*math.exp(offset-(x*x+z*z)/2)/math.sqrt(v)

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", IntegrationWarning)
            value, estimated = quad(f, lo, hi, points=sorted(cuts)[1:-1],
                                    epsabs=1e-11, epsrel=2e-12, limit=500)
        if caught:
            raise RuntimeError(f"Adaptive quad warning {caught[0].message} at {p}")
        result.append(value)
        errors.append(estimated)
    return {"ratio": sum(result), "sides": result, "estimated_errors": errors}


def high_precision_witness(p):
    mp = mpmath.mp
    mp.dps = 65
    m, a, r, eps = map(mp.mpf, (p.m, p.a, p.r, p.epsilon))
    sigma = mp.sqrt(eps)
    target = m+eps
    offset = a*a/(2*target)
    pref = mp.sqrt(target/(2*mp.pi))
    result = []
    for side in (-1, 1):
        shape = mp.mpf(p.shape_on(side))
        cuts = [mp.mpf(0), sigma, 3*sigma, 9*sigma, r]
        if shape > 0:
            for k in range(1, 65):
                cuts.append(r/(1+k/shape))

        def f(u):
            if shape == 0 or u == 0:
                v = eps
            else:
                v = eps+r*r*mp.exp(shape*(1-r/u))
            return pref*mp.exp(offset-((a+side*u)**2+u*u/v)/2)/mp.sqrt(v)

        result.append(mp.quad(f, sorted(set(cuts))))
    return {"ratio": mp.nstr(sum(result), 59), "sides": [mp.nstr(x, 59) for x in result]}


def f32(value):
    return struct.unpack("f", struct.pack("f", value))[0]


def verify_hit(hit):
    raw = hit["point"]
    point = Point(*(f32(raw[k]) for k in ("m","a","r","edge","shape","power","epsilon","delta")),
                  family=raw["family"])
    result = integrate(point)
    margin = point.a*point.a - point.m*math.log(point.m)/(point.m-1)
    return {
        "ratio_gpu": f32(hit["ratio"]),
        "mixture_cpu": result,
        "relative_error": abs(f32(hit["ratio"])/result["ratio"]-1),
        "true_margin": margin,
        "admissible": point.edge < point.m,
    }


def main():
    gpu = json.loads(Path(sys.argv[1]).read_text())
    witnesses = {}
    r = 0.125
    for k, a, family in (("symmetric",1.5,0),("right_plateau",1.25,1)):
        p = Point(2,a,r,r*r,0.25,1,2**-32,0,family)
        witnesses[k] = {
            "input": vars(p), "double": integrate(p), "mpmath_65digits": high_precision_witness(p),
        }
    result = {
        "source": Path(sys.argv[1]).name, "python": sys.version, "scipy": scipy.__version__,
        "mpmath": mpmath.__version__,
        "maxima": {name: verify_hit(row["best"]) for name,row in gpu["summaries"].items()},
        "exact_parameter_witnesses": witnesses,
    }
    Path(sys.argv[2]).write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")


if __name__ == "__main__":
    main()
