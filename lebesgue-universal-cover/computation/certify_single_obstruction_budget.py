"""Evaluate a uniform rational witness size for each fixed area ceiling."""
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

from flint import arb, ctx


def rational(value):
    return arb(value.numerator)/value.denominator


def minimum_budget(ceiling, lower, grid):
    def area(n):
        # Every rounded sampled body has this one-sided support error.
        error = 3*arb.pi()**2/(8*n*n)+arb(2)/(grid+2)
        # A separately sampled disk puts this much inradius into the host.
        disk_error = (1-(arb.pi()/n).cos())/2+arb(2)/(grid+2)
        radius = arb(1)/2-disk_error
        assert radius > 0
        return ceiling*(1+2*error/radius)+arb.pi()*error**2
    lo, hi = 3, 4
    while not area(hi) < lower:
        hi *= 2
    while hi-lo > 1:
        mid = (hi+lo)//2
        if area(mid) < lower:
            hi = mid
        else:
            lo = mid
    assert area(hi) < lower and area(lo) > lower
    return dict(max_points=hi, enlarged_area_upper=str(area(hi)),
                previous_formula_value=str(area(lo)))


if __name__ == "__main__":
    ctx.prec = 160
    grid = 10**10
    raw = Path(__file__).with_name("shallow-cap-audit.json").read_bytes()
    lower = F(json.loads(raw)["area_floor"])
    targets = [".837", ".83730", ".83735", ".83739", ".83740", ".837401", ".83740123"]
    sizes = {t:minimum_budget(rational(F(t)),rational(lower),grid) for t in targets}
    print(json.dumps(dict(grid_denominator=grid+2, area_lower=str(lower),
                          premise_sha256=sha256(raw).hexdigest(), targets=sizes),indent=2))
