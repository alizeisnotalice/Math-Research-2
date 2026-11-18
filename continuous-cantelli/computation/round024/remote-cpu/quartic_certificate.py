"""Exact lower quadrature for a scale-indexed quartic core family."""
import json
import platform
from fractions import Fraction as Q
from math import isqrt
from pathlib import Path

SCALE = 10**9


def sqrt_upper(q):
    units = isqrt(q.numerator*SCALE**2//q.denominator)+1
    return Q(units, SCALE)


def exp_lower(q):
    total, term = Q(1), Q(1)
    for n in range(1, 16):
        term *= -q/n
        total += term
    return max(Q(0), total)


xlo, xhi = Q(6, 5), Q(5, 4)
mlo, mhi = Q(3, 2), Q(2)
klo, khi = Q(15, 4), Q(17, 4)
rhi = Q(1, 8)
xi = Q(1, 1024)
# With c=k/r², both sides satisfy v≤k r²; epsilon=xi r².
assert khi*rhi**2 < mlo
assert xlo > rhi
# log 2 < 7/10 by a positive exponential partial sum.
assert sum(Q(7,10)**i / [1,1,2,6,24][i] for i in range(5)) > 2
# M log M/(M-1) increases for M>1, as M-1-log M>0.
assert xlo**2 > Q(14,10)
# pi>3 by the alternating arctangent series (8 terms).
assert 4*sum(Q((-1)**i, 2*i+1) for i in range(8)) > 3
# The target exponent exceeds 1/3, hence its exponential is <3/4.
assert xlo**2/(2*(mhi+xi*rhi**2)) > Q(1,3)
# The p_{M+epsilon} prefactor is <1/3 since pi>3 and M≥3/2.
assert 6*mlo >= 9
target = Q(1,4)

cells = 1024
accum = 0
for i in range(cells):
    a, b = Q(i,cells), Q(i+1,cells)
    vlo = xi+klo*a**4
    vhi = xi+khi*b**4
    # u=r z, so positive heat kernel scaling cancels with du.
    cost_noise = b*b/(2*vlo)
    cost_source_left = xhi*xhi/2
    cost_source_right = (xhi+rhi*b)**2/2
    density = Q(7,44) * (exp_lower(cost_source_left+cost_noise)
                          +exp_lower(cost_source_right+cost_noise))/sqrt_upper(vhi)
    cell = density/cells
    accum += (SCALE*cell.numerator)//cell.denominator

core_lower = Q(accum,SCALE)
assert core_lower > target
out = {
    "mathematical_ranges": {"kappa":[str(klo),str(khi)], "r":["(0",str(rhi)+"]"],
                            "c":"kappa/r²", "epsilon":"r²/1024", "y":"x",
                            "x":[str(xlo),str(xhi)], "M":[str(mlo),str(mhi)]},
    "rational_cells_per_side":cells,
    "core_density_lower":str(core_lower),
    "target_density_upper":str(target),
    "strict_gap":str(core_lower-target),
    "ratio_lower":str(core_lower/target),
    "python":platform.python_version(),
}
Path('quartic_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
