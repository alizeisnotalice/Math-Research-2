"""Exact rational midpoint-free lower sums for a positive local-density cut."""
from fractions import Fraction as Q
from math import isqrt
import json
import platform
from pathlib import Path

SHRINK = 10**6


def sqrt_lower(q):
    return Q(isqrt(q.numerator*SHRINK**2//q.denominator), SHRINK)


def exp_lower(q):
    # Negative-argument Taylor remainder has the positive sign at odd order.
    total = Q(1)
    term = Q(1)
    for n in range(1, 12):
        term *= -q/n
        total += term
    return max(Q(0), total)


t = Q(1, 3)
t_lower, t_upper = t-Q(1, 10000), t+Q(1, 10000)
def numerator(t):
    return 5913*t**3+954*t**2+96*t-8
def denominator(t):
    return 21*t+6
assert 0 < t_lower < t_upper < Q(2, 3) and numerator(t_lower)>0
# Both numerator and denominator increase on this positive interval.
c2_lower = numerator(t_lower)/denominator(t_upper)
c2_upper = numerator(t_upper)/denominator(t_lower)
c_lower = sqrt_lower(c2_lower)
c_upper = sqrt_lower(c2_upper)+Q(1, SHRINK)
assert c_lower**2 <= c2_lower <= Q(349, 13) <= c2_upper < c_upper**2
d_lower, d_upper = t_lower/c_upper, t_upper/c_lower
x_lower, x_upper = (3*t_lower-2)/c_lower, (3*t_upper-2)/c_upper
eps = Q(1, 2**20)
delta = Q(1, 8)
radius = Q(3, 16)
assert c_upper*radius<1 and d_upper*radius**3<1
# pi<22/7, pi>3 from alternating arctan series.
assert 4*sum(Q((-1)**i, 2*i+1) for i in range(8))>3
k = 256
lower_units = 0
start, end = Q(1, 16), radius
width = (end-start)/k
for i in range(k):
    a, b = start + i*width, start+(i+1)*width
    var_lo = eps + d_lower*a**3
    var_hi = eps + d_upper*b**3
    variance_sqrt_upper = sqrt_lower(var_hi) + Q(1, SHRINK)
    position_abs = max(abs(x_lower+a), abs(x_upper+b))
    difference_abs = max(abs(delta-a), abs(delta-b))
    cost = (position_abs**2 + difference_abs**2/var_lo)/2
    density = Q(7,44)*exp_lower(cost)/variance_sqrt_upper
    lower = width*density
    lower_units += (SHRINK*lower.numerator)//lower.denominator
cubic_lower = Q(lower_units, SHRINK)
# On linear-side u∈[1/16,3/16], the full added variance lies in (3/10,1).
a, b = Q(1,16), radius
assert eps+c_lower*a>Q(3,10) and eps+c_upper*b<1
assert x_lower-b>-Q(2,5) and x_upper-a<0
assert Q(2,5)**2/2+(delta+b)**2/(2*Q(3,10))<Q(1,4)
linear_lower = (b-a)*Q(7,44)*Q(3,4)
lower = cubic_lower+linear_lower
# For every M>1, the target p_{M+eps}(y) is < 1/sqrt(6).
assert 6*lower**2>1
out = dict(
    parameters={"t_interval": [str(t_lower), str(t_upper)], "M": "all M>1",
                "eps": str(eps), "delta":str(delta), "common_radius":str(radius),
                "side_interval": [str(start), str(end)]},
    c_bracket=[str(c_lower),str(c_upper)],
    cubic_cells=k, cubic_density_lower=str(cubic_lower), linear_density_lower=str(linear_lower),
    core_density_lower=str(lower), core_density_lower_decimal=float(lower),
    target_upper_squared="1/6", squared_ratio_lower=str(6*lower**2),
    python=platform.python_version(),
)
Path("exact_uniform.json").write_text(json.dumps(out,indent=2)+"\n")
print(json.dumps(out,indent=2))
