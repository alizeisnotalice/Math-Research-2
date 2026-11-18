"""Exact final arithmetic for the unit mixed-contact local exclusion."""
from fractions import Fraction as Q
from math import factorial


m = 1 + Q(1, 10**6)
r = Q(1, 4)
epsilon = Q(1, 10**12)
delta = Q(1, 10**4)
shift_scale = Q(1, 10**6)
target_variance = m + epsilon
y = 1 + delta
assert 0 < r <= m and r**3 <= m
assert y*y > target_variance

# The normalizing prefactor sqrt((M+epsilon)/(2*pi)) exceeds .39
# by the classical bound pi < 22/7.
normalization = Q(39, 100)
assert normalization**2 * 2 * Q(22, 7) < target_variance

# Left: on 1/100 <= u <= 1/4, (delta+u)^2/(epsilon+u)
# is at most r+2*delta+delta^2/(1/100).
left_a = Q(1, 100)
left_exponent_loss = ((1-left_a)**2 + r + 2*delta + delta**2/left_a - 1)/2
assert 0 < left_exponent_loss < Q(3, 25)
left_scale = Q(99, 100)
assert left_scale**2 * (1 + epsilon/left_a) < 1
left = normalization * (1-Q(3, 25)) * left_scale * 2 * (Q(1, 2) - Q(1, 10))

# Right: on delta-3e-6 <= u <= delta+3e-6, change variables to
# t=(u-delta)/1e-6. Then variance/1e-12 is between 1.9 and 2.1.
right_a = delta - 3*shift_scale
right_b = delta + 3*shift_scale
assert 0 < right_a < right_b < r
assert (epsilon + right_a**3)/shift_scale**2 > Q(19, 10)
assert (epsilon + right_b**3)/shift_scale**2 < Q(21, 10)
right_scale = Q(17, 25)
assert right_scale**2 * Q(21, 10) < 1
assert ((1+right_b)**2 - 1)/2 < Q(1, 1000)

# exp(-s) >= its degree-15 Taylor polynomial for all s >= 0,
# by the Lagrange remainder. Right-end rectangles integrate the
# decreasing function exp(-t^2/(2*1.9)) on [0,3].
def lower_exp(s):
    return sum((-s)**k / factorial(k) for k in range(16))


step = Q(1, 8)
rectangle_sum = 2 * step * sum(lower_exp((j*step)**2 / Q(19, 5)) for j in range(1, 25))
right = normalization * right_scale * Q(999, 1000) * rectangle_sum
assert left + right > Q(11, 10)
print(f"Left lower bound = {left} ({float(left):.15f})")
print(f"Right lower bound = {float(right):.15f}")
print(f"Sum lower bound = {float(left+right):.15f} > 11/10 (exact comparison)")
