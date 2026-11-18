"""Rational verification of the sign of the critical first correction."""
from fractions import Fraction as Q
from math import factorial


def sin_polynomial(t, degree):
    return sum((-1)**k * t**(2*k+1)/factorial(2*k+1) for k in range((degree+1)//2))


# For 0<x<pi/2, use monotonicity and rational pi bounds.
pi_lo, pi_hi = Q(314159, 100000), Q(22, 7)
def atan_partial(x, nterms):
    return sum((-1)**k*x**(2*k+1)/(2*k+1) for k in range(nterms))


# Machin's identity and alternating rational series.
q, w = Q(1, 5), Q(1, 239)
assert 4*(4*atan_partial(q, 8)-atan_partial(w, 9)) > pi_lo
assert 4*(4*atan_partial(q, 9)-atan_partial(w, 8)) < pi_hi
for n in (4,6,8):
    lo, hi = pi_lo/(2*(n-1)), pi_hi/(2*(n-1))
    slo, shi = sin_polynomial(lo, 7), sin_polynomial(hi, 5)
    tlo, thi = sin_polynomial(2*lo, 7), sin_polynomial(2*hi, 5)
    ell, upper = Q(n-2, 2)*slo, Q(n-2, 2)*shi
    if n < 8:
        assert tlo**n > upper**(n+1)
        sign = "positive"
    else:
        assert thi**n < ell**(n+1)
        sign = "negative"
    print(n, sign, float(slo), float(shi), float(tlo), float(thi))
