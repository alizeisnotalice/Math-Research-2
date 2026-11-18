"""Independent optimization and moving-saddle checks for the Metal profiles."""
import cmath
import json
import math
from pathlib import Path
import struct
import sys

import mpmath as mp
import scipy
from scipy.optimize import brentq


def f32(value):
    return struct.unpack("f", struct.pack("f", value))[0]


def unpack(inp):
    return int(inp["n"]), *(f32(inp[k]) for k in ("logEta", "r", "hR", "bR2", "a", "logLambda", "logK"))


def numeric_profile(n, logeta, r, h_r, b_r2, a, loglam, logk):
    s = (2.0**logeta)*r
    lam = 2.0**loglam
    h, b = h_r/r, b_r2/r**2
    theta = math.pi/(2*(n-1))
    A = (1-1/n)*math.sin(theta)
    z0 = (((n-2)/2)*math.sin(theta)*math.exp(n*logk/(2*n-2)))**(1/n)
    kappa = z0**(2*n-2)/(n-2)
    L = r/s

    def minimize(sign, source_weight):
        H = sign*a*s/lam if source_weight else 0.0
        J = s*s/lam if source_weight else 0.0
        u, v = sign*h*s, b*s*s

        def derivative(z):
            return (-(n-2)*kappa/z**(n-1)
                    + z**(n-1)*(1+((n+1)/n)*u*z+((n+2)/n)*v*z*z)+H+J*z)

        loc = brentq(derivative, min(z0/8, L/2), L, xtol=1e-14) if derivative(L) > 0 else L
        act = kappa/loc**(n-2) + (loc**n/n)*(1+u*loc+v*loc**2) + H*loc+J*loc**2/2
        return act, loc

    act, loc = min(minimize(-1, False), minimize(1, False))
    source, _ = min(minimize(-1, True), minimize(1, True))
    first_real = -(abs(h)*s/n)*z0**(n+1)
    first_analytic = -(abs(h)*s/n)*math.sin(2*theta)
    second_analytic = -(s*s)*math.sin(3*theta)*(b/n - h*h*(n+1)**2/(2*n*n*(n-1)))
    base = A*math.exp(n*logk/(2*n-2))

    def moving_decay(z):
        u, v = h*s, b*s*s
        for step in range(1, 17):
            fraction = step/16
            U, V = u*fraction, v*fraction*fraction
            for _ in range(12):
                f = z**(n-1)*(1+(n+1)*U*z/n+(n+2)*V*z*z/n)-1j
                fp = z**(n-2)*((n-1)+(n+1)*U*z+(n+2)*(n+1)*V*z*z/n)
                dz = f/fp
                z -= dz
                if abs(dz) < 1e-14:
                    break
        G = z**n*(1+u*z+v*z*z)/n - 1j*z
        residual = z**(n-1)*(1+(n+1)*u*z/n+(n+2)*v*z*z/n)-1j
        return G.real, abs(residual)

    z = cmath.exp(1j*theta)
    saddles = [moving_decay(v) for v in (z, -z.conjugate())]
    analytic = min(value for value, error in saddles)
    return {
        "t": loc/z0-1,
        "realGap": act-A,
        "firstGap": base-A+first_real-first_analytic,
        "analyticGap": act-A-first_analytic-second_analytic,
        "sourceGap": source-A,
        "remainder": (act-base-first_real)/z0**n,
        "logAmplitude": math.log(act)+math.log(lam),
        "exact_moving_saddle_gap": act-analytic,
        "max_moving_saddle_residual": max(error for value, error in saddles),
        "R0": base,
        "A": A,
    }


def high_precision(n, logeta, r, h_r, b_r2, a, loglam, logk):
    mp.mp.dps = 65
    logeta, r, h_r, b_r2, a, loglam, logk = map(mp.mpf, (logeta, r, h_r, b_r2, a, loglam, logk))
    s, lam = mp.power(2, logeta)*r, mp.power(2, loglam)
    h, b = h_r/r, b_r2/r**2
    theta = mp.pi/(2*(n-1))
    A = (1-mp.mpf(1)/n)*mp.sin(theta)
    z0 = (((n-2)/mp.mpf(2))*mp.sin(theta)*mp.exp(n*logk/(2*n-2)))**(mp.mpf(1)/n)
    kappa = z0**(2*n-2)/(n-2)

    def act(z, sign, weight):
        out = kappa/z**(n-2)+z**n*(1+sign*h*s*z+b*s*s*z*z)/n
        if weight:
            out += (sign*a*s*z+s*s*z*z/2)/lam
        return out

    def minimum(sign, weight):
        derivative = lambda z: mp.diff(lambda zz: act(zz, sign, weight), z)
        L = r/s
        if derivative(L) <= 0:
            return act(L, sign, weight)
        root = mp.findroot(derivative, (z0*mp.mpf(".97"), z0*mp.mpf("1.03")))
        return act(root, sign, weight)

    actual = min(minimum(-1, False), minimum(1, False))
    source = min(minimum(-1, True), minimum(1, True))
    first_analytic = -(abs(h)*s/n)*mp.sin(2*theta)
    second_analytic = -(s*s)*mp.sin(3*theta)*(b/n - h*h*(n+1)**2/(2*n*n*(n-1)))
    return {
        "realGap": mp.nstr(actual-A, 50),
        "analyticGap": mp.nstr(actual-A-first_analytic-second_analytic, 50),
        "sourceGap": mp.nstr(source-A, 50),
    }


def main():
    gpu = json.loads(Path(sys.argv[1]).read_text())
    comparisons = {}
    for name,row in gpu["summaries"].items():
        if row["best"] is None:
            continue
        raw = row["best"]
        args = unpack(raw["input"])
        actual = numeric_profile(*args)
        precise = high_precision(*args)
        errors = {
            k: abs(f32(raw["answer"][k])-float(precise.get(k, actual[k])))
            for k in ("t","realGap","firstGap","analyticGap","sourceGap","remainder","logAmplitude")
        }
        comparisons[name] = {
            "input64": args,
            "metal": raw["answer"],
            "double": actual,
            "mpmath_65digits": precise,
            "absolute_error": errors,
        }
    canon = {}
    for n in (4, 6, 8):
        for h_r in (0.0, .2):
            for eta in (-4.0, -9.0, -14.0):
                args = (n, eta, 0.25, h_r, .1, 2.0, 8.0 + (n-2)*(-eta), 0.0)
                canon[f"n{n}_hR{h_r}_logeta{eta}"] = numeric_profile(*args)
    result = {
        "source": Path(sys.argv[1]).name, "python": sys.version,
        "scipy": scipy.__version__, "mpmath": mp.__version__,
        "comparison": comparisons, "canonical": canon,
    }
    Path(sys.argv[2]).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
