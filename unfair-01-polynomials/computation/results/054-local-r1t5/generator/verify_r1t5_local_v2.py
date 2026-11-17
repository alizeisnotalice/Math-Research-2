#!/Users/zzh/.venvs/main/bin/python
"""Fresh exact certificate for an algebraic-integer r=1,t=5 local prefix.

It deliberately does not modify the earlier sqrt(3)/2 certificate.  The
assertions here concern only C_0,...,C_5; C_6 is exhibited as non-binary.
"""
from fractions import Fraction
import hashlib
import json
import os
import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "r1t5-local-certificate-v2.json"
z = sp.Symbol("z")


def conv(A, B):
    C = [0] * (len(A) + len(B) - 1)
    for i, ai in enumerate(A):
        for j, bj in enumerate(B):
            C[i + j] += ai * bj
    return [sp.cancel(x) for x in C]


def family(u):
    v = 1 - u
    s = 1 - u * v
    t = 1 - u * s
    rr = 1 - v * t
    b = sp.cancel((1 - v * rr - t * s) / (2 * u - 1))
    return (
        [1, u, 0, t, sp.cancel(rr - b), 0],
        [1, v, s, 0, b, 0],
        {"u": u, "v": v, "s": s, "t": t, "r": rr, "b": b, "A4": sp.cancel(rr-b)},
    )


def qstr(x):
    return str(sp.expand(sp.radsimp(sp.cancel(x))))


def rat(x):
    x = sp.Rational(x)
    return [int(x.p), int(x.q)]


def linear_pair(x, radical):
    q = sp.Poly(sp.expand(sp.radsimp(sp.cancel(x))), radical)
    assert q.degree() <= 1, (x, q)
    return q.nth(0), q.nth(1)


def interval_bounds(x, lo, hi, radical):
    """Exact lower/upper bounds for p+q*radical on (lo,hi)."""
    p, q = linear_pair(x, radical)
    a, b = p + q * lo, p + q * hi
    return min(a, b), max(a, b)


def minpoly_string(x):
    return str(sp.Poly(sp.minpoly(x, z), z).as_expr())


def main():
    radical = sp.sqrt(5)
    u = (-5 + 3 * radical) / 2
    A, B, p = family(u)
    C = conv(A, B)

    assert all(sp.simplify(C[i] - 1) == 0 for i in range(6))
    c6 = C[6]
    assert sp.simplify(c6 - p["A4"] * p["s"]) == 0
    assert sp.simplify(c6 - (-31770 + 14208 * radical)) == 0

    # A tight purely rational isolating interval, certified by endpoint squares.
    lo, hi = sp.Rational(22360679, 10**7), sp.Rational(22360680, 10**7)
    assert lo > 0 and lo**2 < 5 < hi**2

    named = {
        "u": p["u"], "v": p["v"], "s": p["s"], "t": p["t"],
        "b": p["b"], "A4": p["A4"], "C6": c6,
    }
    expected_minpoly = {
        "u": "z**2 + 5*z - 5",
        "v": "z**2 - 7*z + 1",
        "s": "z**2 - 42*z + 36",
        "t": "z**2 - 242*z + 61",
        "b": "z**2 + 110*z - 100",
        "A4": "z**2 + 1545*z - 95",
        "C6": "z**2 + 63540*z - 3420",
    }
    positivity = {}
    for name, value in named.items():
        lower, upper = interval_bounds(value, lo, hi, radical)
        assert lower > 0 and upper < 1, (name, lower, upper)
        mp = minpoly_string(value)
        assert mp == expected_minpoly[name], (name, mp)
        positivity[name] = {
            "linear_in_sqrt5": [qstr(x) for x in linear_pair(value, radical)],
            "strict_lower_bound": rat(lower),
            "strict_upper_bound": rat(upper),
            "minimal_polynomial": mp,
        }

    # The two zero slots in A/B are intentional.  Every nonzero, nonunit
    # coefficient is in the preceding certified (0,1) list.
    obj = {
        "status": "passed_exact_algebraic_integer_local_prefix_only",
        "scope": (
            "Exact local equations C_0,...,C_5=1 for a necessary r=1,t=5 "
            "prefix. This is neither a full nonnegative factorization nor an "
            "original-problem counterexample, and it does not establish global "
            "r=1,t=5 feasibility."
        ),
        "u": qstr(u),
        "sqrt5_isolating_interval": [rat(lo), rat(hi)],
        "endpoint_square_certificate": {
            "lo_squared_lt_5": rat(5 - lo**2),
            "hi_squared_minus_5": rat(hi**2 - 5),
        },
        "A": [qstr(x) for x in A],
        "B": [qstr(x) for x in B],
        "parameters": {name: qstr(value) for name, value in p.items()},
        "C_0_to_6": [qstr(x) for x in C[:7]],
        "nonzero_nonunit_coefficients_and_C6": positivity,
        "first_nonbinary_coefficient": "C_6",
        "C6": qstr(c6),
        "software": {"python": sys.version, "sympy": sp.__version__, "interpreter": sys.executable},
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    tmp = OUT.with_suffix(".partial")
    with tmp.open("w") as out:
        json.dump(obj, out, sort_keys=True, indent=2)
        out.write("\n")
        out.flush()
        os.fsync(out.fileno())
    os.replace(tmp, OUT)
    print(json.dumps({"status": obj["status"], "C6": qstr(c6), "sympy": sp.__version__}, sort_keys=True))


if __name__ == "__main__":
    main()
