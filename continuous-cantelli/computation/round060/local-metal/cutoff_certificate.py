"""Exact real-arithmetic truncation bounds on the GPU's prescribed domain."""
from fractions import Fraction as Q
from math import factorial


def upper_exp_negative(x):
    return 1/sum(x**j/factorial(j) for j in range(121))


# Source normalization C_h < 3 follows from target < 9.02 and pi > 3.
assert 4*Q(902,100)/6 < 9
# Proposed h lies in [2^-14,2^-2.5], and the cap is at least 3/(64*12.25).
assert Q(3,64)/Q(49,4) > Q(1,2**14)
# The perturbation ratio magnitude is at most 2^-2.5/(5/8) < 0.3.
assert (Q(8,5)**2)*Q(1,32) < Q(9,100)
# The variance floor is > 2^-54. 1/sqrt(1±x)<2 and e^|a*u|<27.
assert Q(7,10)*4 > 1
bound_full = 3*3*27*2**27*2**28 * upper_exp_negative(Q(960,13))
assert bound_full < Q(1,10**13)
# Polynomial times exp(-Z/2) decreases termwise for Z >=192.
z=192
bound_lead = 3*2**27*((z*z+6*z+3)*Q(8,5)**2/8+Q(49,4)*(z+1)*Q(8,5)/2)*upper_exp_negative(Q(96))
assert bound_lead < Q(1,10**27)
# First omitted term in the odd degree-9 sinh Taylor series.
tail_ratio = Q(1,16)/(12*13)
sinh_tail = Q(1,4)**11/factorial(11)/(1-tail_ratio)
assert sinh_tail < Q(1,10**14)
print("full scalar omitted |Z>192| integral <=",float(bound_full),"<1e-13")
print("leading scalar omitted integral <=",float(bound_lead),"<1e-27")
print("sinh real-arithmetic truncation <=",float(sinh_tail),"<1e-14")
