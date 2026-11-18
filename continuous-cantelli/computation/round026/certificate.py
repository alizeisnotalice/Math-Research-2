"""Fraction checks for the finite two-sided tilted plateau witness."""
from fractions import Fraction as Q
from math import factorial
import json
import platform
from pathlib import Path

lam, a, m, b, r = Q(20), Q(1), Q(3,4), Q(2), Q(1)
wide, narrow = Q(49,100), Q(97,200)
assert (1-wide)*lam>r
assert a*r*r>=m
# e<11/4 by bounding terms n>=5 with a ratio-1/6 geometric tail.
assert sum(Q(1,factorial(i)) for i in range(5))+Q(1,100)<Q(11,4)
# sqrt(2pi) in (12/5,8/3) by 3<pi<22/7.
assert 4*sum(Q((-1)**i,2*i+1) for i in range(8))>3
assert Q(12,5)**2<6 and Q(8,3)**2>Q(44,7)
# On |x|≤r, bound the core by the all-real quadratic integral.
# exp(1/(2a))/sqrt(1+a lambda²) < 2/lambda.
assert a==1 and Q(11,4)<4
core_upper = 2/lam
# On the high plateau use e^{-(b-1)lambda²/2} < 2^-200 <1/1000.
assert (b-1)*lam*lam/2==200 and Q(1,2**200)<Q(1,1000)
plateau_upper = Q(1,1000)
# Gaussian upper Mills bound, with positive exponent action <2.
q = lam*wide
action_wide = (1-m-wide*wide)*lam*lam/2
assert 0<action_wide<2
outer_upper = 2*Q(11,4)**2/(q*Q(12,5))
upper = core_upper+plateau_upper+outer_upper
assert upper < 1
# Conversely the saturated discontinuous envelope at the narrower window
# has F>1. Keep both outer tails and subtract an upper bound for the
# doubly counted tilted core mass.
q = lam*narrow
action_narrow = (1-m-narrow*narrow)*lam*lam/2
exp_action_lower = sum(action_narrow**i/factorial(i) for i in range(6))
outer_lower = 2*exp_action_lower*(q/(q*q+1))/Q(8,3)
core_action = (1-m)*lam*lam/2 - (lam-r)**2/2
assert core_action < -10
lower = outer_lower-Q(1,2**10)
assert lower > 1
out = dict(
    lambda_=str(lam), quadratic_lower_a=str(a), core_radius=str(r),
    off_core_lower_m=str(m), high_plateau_lower_b=str(b),
    symmetric_distance_wider=str(wide), upper_wider=str(upper),
    symmetric_distance_narrower=str(narrow), lower_for_narrower_saturated=str(lower),
    strict_gaps=[str(1-upper),str(lower-1)], python=platform.python_version())
Path('certificate.json').write_text(json.dumps(out,indent=2)+"\n")
print(json.dumps(out,indent=2))
