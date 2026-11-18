"""Fraction checks for the right-core density obstruction suggested by Metal.

For a candidate exact kernel let c=sqrt(C), d=t/c, x=(3*t-2)/c.
Take the assumed right core v(x+u)=d*u**3 for 0<u<r. At
 y=x-a/c and smoothing variance e, put z=c*u. With Q=3*t-2,
 V(z)=e+t*z**3/C**2, the right-core density divided by p_(M+e)(y)
 per unit z is sqrt((M+e)/(2*pi*C*V(z))) * exp(E(z)), where
 E(z)=(Q-a)**2/(2*C*(M+e))-(Q+z)**2/(2*C)-(z+a)**2/(2*C*V(z)).
On l<=z<=h<-Q the input square is bounded by (Q+l)**2, the
 displacement square by (h+a)**2, and V by its endpoint values.
The tests below verify the rational bounds. Positivity of the rest of
 the mixture makes the certificate independent of its contribution.
"""
import json
from fractions import Fraction as F
from math import factorial


t = F(50421, 10**6)
c_squared = (5913 * t**3 + 954 * t**2 + 96 * t - 8) / (21 * t + 6)
q = 3 * t - 2
radius = F(1, 4)
m = F(10001, 10000)
e = F(1, 100000)
a = F(1, 125)
left = F(139, 10000)
right = F(144, 10000)


def variance(z):
    return e + t * z**3 / c_squared**2


assert c_squared > 0 and 0 < left < right < -q
assert right**2 < c_squared * radius**2
assert c_squared * radius**2 < m**2
assert t**2 * radius**6 < c_squared * m**2
lower_exponent = ((q-a)**2 / (2*c_squared*(m+e))
                  - (q+left)**2 / (2*c_squared)
                  - (right+a)**2 / (2*c_squared*variance(left)))
assert lower_exponent > 5
lower_prefactor_squared = (right-left)**2 * (m+e) / (8*c_squared*variance(right))
assert lower_prefactor_squared > F(1, 40)**2
lower_ratio = sum(F(5**j, factorial(j)) for j in range(5)) / 40
assert lower_ratio == F(523, 320) > 1
print(json.dumps({
    "t": str(t), "c_squared": str(c_squared), "radius": str(radius),
    "M": str(m), "epsilon": str(e), "a": str(a),
    "z_left": str(left), "z_right": str(right),
    "right_core_ratio_strictly_exceeds": str(lower_ratio)
}, indent=2))
