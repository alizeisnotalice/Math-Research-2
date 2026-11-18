"""Exact finite quartic example outside the positional forbidden interval.

Let v(a+u)=c*u**4 for |u|<r and assume 0<=v<=M globally.
With epsilon smoothing and observer y, evaluate the left input z=a-w
for l<=w<=h<r. Write V(w)=epsilon+c*w**4. Its relative density per unit
w is sqrt((M+epsilon)/(2*pi*V(w))) times exp(E(w)), where
E(w)=y*y/(2*(M+epsilon))-(a-w)**2/2-(y-a+w)**2/(2*V(w)).
The endpoint bounds are exact and the rest of the mixture is nonnegative.
"""
import json
from fractions import Fraction as F
from math import factorial

a = F(5)
c = F(3, 16)
r = F(3, 2)
m = F(1025, 1024)
epsilon = F(1, 1024)
y = F(55, 8)
l = F(23, 16)
h = F(95, 64)


def variance(w):
    return epsilon + c * w**4


assert 0 < l < h < r < a < y
assert c * r**4 < 1 < m < a**2
# Gaussian core mass <= Gaussian tail at a-r. By Mills' bound, pi>1/2,
# (a-r)**2/2>6 and exp(1)>8/3, this tail is less than tail_upper.
assert (a - r)**2 / 2 > 6
assert sum(F(1, factorial(j)) for j in range(4)) == F(8, 3)
tail_upper = F(3, 8)**6 / (a - r)
assert tail_upper < m - 1 and m * tail_upper < 1
# a**2 > M implies a**2 > M*log(M)/(M-1) since log(M)<M-1.
lower_exponent = y**2 / (2*(m+epsilon)) - (a-l)**2 / 2 - (y-a+h)**2 / (2*variance(l))
assert lower_exponent > 10
# Use pi<4, increasing V, and positive prefactor.
lower_prefactor_squared = (h-l)**2 * (m+epsilon) / (8*variance(h))
assert lower_prefactor_squared > F(1, 64)**2
lower_ratio = sum(F(10**j, factorial(j)) for j in range(5)) / 64
assert lower_ratio == F(1933, 192) > 1
print(json.dumps({
    "a": str(a), "c": str(c), "r": str(r), "M": str(m),
    "epsilon": str(epsilon), "y": str(y), "left_abs_u_lower": str(l),
    "left_abs_u_upper": str(h), "core_mass_upper": str(tail_upper),
    "core_ratio_strictly_exceeds": str(lower_ratio)
}, indent=2))
