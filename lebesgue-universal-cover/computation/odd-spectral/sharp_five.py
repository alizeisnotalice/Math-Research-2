"""Enclose the sharp first rotation moment when odd frequencies above five vanish."""
from fractions import Fraction as F
import json
from pathlib import Path
import sys
from flint import arb, ctx

def polynomial(x):
    return 64*x**6 - 240*x**4 - 40*x**3 + 180*x**2 + 48*x - 5

def interval(x):
    return arb(x.numerator) / x.denominator

def run():
    ctx.prec = 256
    lo, hi = F("0.0805"), F("0.0806")
    assert polynomial(lo) < 0 < polynomial(hi)
    derivative_lower = 48 + 360*lo - 120*hi**2 - 960*hi**3
    assert derivative_lower > 0
    for _ in range(170):
        middle = (lo + hi) / 2
        if polynomial(middle) < 0:
            lo = middle
        else:
            hi = middle
    assert polynomial(lo) < 0 < polynomial(hi)
    d = interval(lo).union(interval(hi))
    u = ((1 + 6*d - 2*d**3)/(24*d)).sqrt()
    a, b = u - d/2, u + d/2
    assert 0 < a < b < 1
    c = arb.pi() / (4*(1 - 2*d))
    assert arb("0.9361264544688870") < c < arb("0.9361264544688872")
    # The last identity is an integer polynomial identity after eliminating u^2.
    # Check its two sides as polynomials using rational coefficient arithmetic.
    from flint import fmpq_poly
    x = fmpq_poly([0, 1])
    numerator = 1 + 6*x - 2*x**3
    # 16 (5d u^4 + 5 d^3 u^2/2 + d^5/16) - (3+10d), times 36d.
    eliminated = 5*numerator**2 + 60*x**3*numerator + 36*x**6 - 36*x*(3+10*x)
    assert eliminated == -polynomial(x)
    return dict(precision=ctx.prec, root_rational_interval=[str(lo), str(hi)],
                first_change=str(a), second_change=str(b), first_rotation_moment=str(c),
                area_factor=str(c*c))

if __name__ == "__main__":
    data = json.dumps(run(), indent=2) + "\n"
    if len(sys.argv) == 2:
        Path(sys.argv[1]).write_text(data)
    else:
        print(data, end="")
