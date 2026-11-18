"""Independent double and multiprecision checks for quartic Metal output."""
import json
import math
from pathlib import Path
import struct
import sys

import mpmath
import scipy
from scipy.integrate import quad


def panels(r, c, eps, delta):
    points = {0.0, r}
    for scale in (math.sqrt(eps), 1/math.sqrt(c), eps**.25/c**.25, abs(delta)):
        for factor in (1/64, 1/16, 1/4, 1/2, 1, 2, 4, 8, 16, 64):
            t = factor*scale
            if 0 < t < r:
                points.add(t)
    if delta > 0:
        scale = math.sqrt(eps+c*delta**4)
        for k in (-16, -8, -4, -2, -1, 0, 1, 2, 4, 8, 16):
            t = delta+k*scale
            if 0 < t < r:
                points.add(t)
    return sorted(points)


def integrate(m, a, c, r, eps, delta):
    target = m + eps
    offset = (a + delta)**2/(2*target)
    prefactor = math.sqrt(target/(2*math.pi))
    sides = []
    for sign in (-1, 1):
        def kernel(u):
            x = a+sign*u
            shift = delta-sign*u
            v = eps+c*u**4
            return prefactor/math.sqrt(v)*math.exp(offset-(x*x+shift*shift/v)/2)

        sides.append(quad(kernel, 0, r, points=panels(r, c, eps, sign*delta)[1:-1],
                          epsabs=1e-11, epsrel=2e-12, limit=250))
    return {
        "ratio": sum(v for v, err in sides),
        "errors": [err for v, err in sides],
        "sides": [v for v, err in sides],
    }


def f32(value):
    return struct.unpack("f", struct.pack("f", value))[0]


def verify(hit):
    args = [f32(hit[key]) for key in ("m", "a", "c", "r", "epsilon", "delta")]
    m, a, c, r, eps, delta = args
    result = integrate(*args)
    return {
        "gpu": hit["ratio"],
        "cpu": result,
        "relative_error": abs(f32(hit["ratio"])/result["ratio"]-1),
        "true_margin": a*a-m*math.log(m)/(m-1),
        "admissible": c*r**4 <= m and m > 1,
    }


def high_precision():
    mp = mpmath.mp
    mp.dps = 65
    m, a, c, r, eps = mp.mpf(2), mp.mpf(3)/2, mp.mpf(4096), mp.mpf(1)/8, mp.mpf(1)/2**32
    sigma = mp.sqrt(eps)
    pref = mp.sqrt((m+eps)/(2*mp.pi))

    def kernel(u):
        v = eps+c*u**4
        return (2*pref/mp.sqrt(v)*mp.exp(a*a/(2*(m+eps))-(a*a+u*u+u*u/v)/2)
                *mp.cosh(a*u))

    points = [0, sigma, 3*sigma, 16*sigma, 1/mp.sqrt(c), r]
    return mp.nstr(mp.quad(kernel, points), 60)


def main():
    gpu = json.loads(Path(sys.argv[1]).read_text())
    report = {
        "python": sys.version,
        "scipy": scipy.__version__,
        "mpmath": mpmath.__version__,
        "input": Path(sys.argv[1]).name,
        "maxima": {name: verify(value["best"]) for name, value in gpu["summaries"].items()},
        "rational_input_check": {
            "float64": integrate(2, 1.5, 4096, .125, 2**-32, 0),
            "mpmath_65digits": high_precision(),
        },
    }
    Path(sys.argv[2]).write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
