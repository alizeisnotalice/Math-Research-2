"""Rational arithmetic for an off-center symmetric quartic density exclusion."""
from fractions import Fraction as Q
from math import factorial


a = Q(3, 2)
m = Q(2)
r = Q(1, 8)
c = Q(4096)
eps = Q(1, 2**32)
sigma = Q(1, 2**16)
root_c = Q(64)

# Admissibility and strict escape from the previous central exclusion.
# The old threshold in a^2 is 2*log(2) < 2.
assert a*a > 2
assert sigma**2 == eps and root_c**2 == c
assert c*r**4 < m
assert 3*sigma < 1/root_c < r


def exp_lower(s, odd_degree):
    return sum((-s)**k / factorial(k) for k in range(odd_degree+1))


variance = m + eps
normalization = Q(14, 25)
assert normalization**2 * 2 * Q(22, 7) < variance   # π < 22/7
# For x=a±u, summing the source factors yields the factor cosh(a*u)≥1.
# Thus the common exponent is bounded below by -4/7 on [0,r].
loss = a*a*(1-1/variance)/2 + r*r/2
assert loss < Q(4, 7)
assert exp_lower(Q(4, 7), 3) > normalization

relative_scale = Q(999, 1000)
assert relative_scale**2 * (1 + c * (3*sigma)**4/eps) < 1
assert relative_scale**2 * (1 + eps/(c*(1/root_c)**4)) < 1

# Near window: u=σt, 0≤t≤3. Using variance at least ε yields a
# Gaussian exponential exp(-t²/2) as the lower bound. Integrate
# via right-end rectangles and a pointwise Taylor lower bound.
step = Q(1, 8)
near = 2*relative_scale*step*sum(
    exp_lower((j*step)**2/2, 31) for j in range(1, 25)
)
# Broad window: 1/64≤u≤r, t=1/(64u). Using variance at least
# c u⁴ yields an integrand ≥(999/1000) exp(-t²/2) after substitution.
def broad_antiderivative(t):
    return t-t**3/6


broad = 2*relative_scale*(broad_antiderivative(Q(1)) -
                          broad_antiderivative(1/(root_c*r)))
lower = normalization**2 * (near+broad)
assert lower > Q(7, 6)
print(f"Near contribution to normalized integral lower bound: {float(near):.15f}")
print(f"Broad contribution to normalized integral lower bound: {float(broad):.15f}")
print(f"Core/target density ratio lower bound: {float(lower):.15f} > 7/6 (exact comparison)")
