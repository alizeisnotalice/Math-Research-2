"""Fraction checks for discarded heat integrals and stable leakage calculations."""
from fractions import Fraction as Q
from math import factorial


def exp_lower(x, order):
    return sum(x**k/factorial(k) for k in range(order+1))


q_min = Q(7,20)
tail_rate = 37/q_min
exp37 = exp_lower(Q(37), 100)
raw_tail = 1/exp37
leakage_tail = tail_rate/exp37
assert raw_tail < Q(1,10**16)
assert leakage_tail < Q(1,10**14)

response_error_per_eta = 2/(q_min*q_min*exp_lower(Q(2), 18))
assert response_error_per_eta < Q(9,4)

# For |shift| < .1, the subtraction branch uses exp(shift) - 1's
# degree-6 polynomial. Absorb one shift into rate*eta and use
# max rate*exp(-.6*rate) = 1/(.6*e).
expminus_scaled_tail = Q(10,9)*Q(1,10**6)/factorial(7)/(Q(3,5)*Q(8,3))
assert expminus_scaled_tail < Q(14,10**11)
print("Raw coefficient omitted tail <",float(raw_tail))
print("Scaled leakage and derivative omitted tail <",float(leakage_tail))
print("Linear-response remainder / eta <",float(response_error_per_eta))
print("Small-shift normalized expm1 polynomial error <",float(expminus_scaled_tail))
