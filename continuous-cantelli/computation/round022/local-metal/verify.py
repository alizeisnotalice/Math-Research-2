import json
from math import exp, isfinite, pi, sqrt
from pathlib import Path
import sys

import mpmath
import scipy
from scipy.integrate import quad


def breakpoints(r, eps, delta, right):
    points = {0.0, r}
    for p in range(-15, 2):
        t = 10.0**p
        if 0 < t < r:
            points.add(t)
    if right and delta > 0:
        s = sqrt(eps + delta**3)
        for k in (-16, -8, -4, -2, -1, 0, 1, 2, 4, 8, 16):
            t = delta + k*s
            if 0 < t < r:
                points.add(t)
    return sorted(points)


def integrate(m, r, eps, delta):
    target_variance = m + eps
    offset = (1.0 + delta)**2 / (2.0 * target_variance)
    prefactor = sqrt(target_variance / (2.0 * pi))
    pieces = []
    for right in (False, True):
        def term(u):
            x = 1 + u if right else 1 - u
            z = delta - u if right else delta + u
            variance = eps + (u**3 if right else u)
            return prefactor / sqrt(variance) * exp(offset - (x*x + z*z / variance)/2)

        points = breakpoints(r, eps, delta, right)
        pieces.append(quad(term, 0, r, points=points[1:-1], epsabs=2e-11, epsrel=2e-12, limit=300))
    return {
        "ratio": sum(piece[0] for piece in pieces),
        "side_values": [piece[0] for piece in pieces],
        "estimated_error": sum(piece[1] for piece in pieces),
    }


def high_precision_witness():
    mp = mpmath.mp
    mp.dps = 65
    m, r, eps, delta = map(mp.mpf, ("1.000001", ".25", ".000000000001", ".0001"))
    v = m + eps
    prefactor = mp.sqrt(v / (2*mp.pi))
    pieces = []
    for right in (False, True):
        def term(u):
            x = 1 + u if right else 1 - u
            z = delta - u if right else delta + u
            variance = eps + (u**3 if right else u)
            return prefactor / mp.sqrt(variance) * mp.exp(((1 + delta)**2/v - x*x - z*z/variance)/2)

        points = [mp.mpf(0), *(mp.mpf(10)**p for p in range(-15, 0)), r]
        if right:
            scale = mp.sqrt(eps + delta**3)
            points.extend(delta + k*scale for k in (-16, -8, -4, -2, 0, 2, 4, 8, 16))
        pieces.append(mp.quad(term, sorted(t for t in set(points) if 0 <= t <= r)))
    return {
        "ratio": mp.nstr(sum(pieces), 60),
        "side_values": [mp.nstr(piece, 60) for piece in pieces],
        "decimal_precision": mp.dps,
    }


def verify_hit(hit):
    m, r, eps, d = (hit[key] for key in ("M", "r", "epsilon", "delta"))
    # Swift's JSON emits short round-trip binary32 decimals. Reconstitute the
    # device inputs to avoid treating a rounded displayed radius as exact.
    import struct
    f32 = lambda value: struct.unpack("f", struct.pack("f", value))[0]
    m, r, eps, d = map(f32, (m, r, eps, d))
    actual = integrate(m, r, eps, d)
    assert isfinite(actual["ratio"])
    return {
        "input": {"M": m, "r": r, "epsilon": eps, "delta": d},
        "admissible_radius": r <= min(m, m**(1.0 / 3.0)),
        "gpu_ratio": f32(hit["ratio"]),
        "verified": actual,
        "relative_error": abs(f32(hit["ratio"])/actual["ratio"] - 1),
    }


def main():
    gpu = json.loads(Path(sys.argv[1]).read_text())
    report = {
        "python": sys.version,
        "scipy": scipy.__version__,
        "mpmath": mpmath.__version__,
        "gpu_input": Path(sys.argv[1]).name,
        "verified_gpu_maxima": {
            name: verify_hit(summary["hit"]) for name, summary in gpu["summaries"].items()
            if summary["hit"] is not None
        },
        "exact_decimal_witness": {
            "input": {"M": 1.000001, "r": 0.25, "epsilon": 1e-12, "delta": 0.0001},
            "double": integrate(1.000001, 0.25, 1e-12, 0.0001),
            "high_precision": high_precision_witness(),
        },
    }
    Path(sys.argv[2]).write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
