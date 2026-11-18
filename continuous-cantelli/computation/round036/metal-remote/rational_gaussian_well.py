"""Uniform integrated-transform certificate for a shallow stretched Gaussian well.

For any a>=16, assume  f(a+y) <= 1-(11/10)/a**2*exp(-4*y**2)
for |y|<=3/2. In J=a+[-1/8,1/8], the original action is uniformly at
least (11/20)(127/128)^2*exp(-4*y**2) - (|y|+1/8)^2/2.
The two terms decrease in |y|. Rational right-endpoint bounds and Taylor
lower sums evaluate the corresponding rectangle lower bound.
"""
from fractions import Fraction as F
import json
from math import factorial

MINIMUM_CENTER = F(16)
DEPTH = F(11, 10)
CORE = F(3, 2)
HALF_WIDTH = F(1, 8)
PANELS = 128
SCALE = 2**24
NORMAL_LOWER = F(199, 500)


def exp_lower(argument, degree):
    # Odd degree leaves a nonnegative Lagrange remainder for all real arguments.
    return sum(argument**j/factorial(j) for j in range(degree+1))


def main():
    assert NORMAL_LOWER**2 * F(44, 7) < 1  # pi < 22/7
    minimum_scaled_t = 1-HALF_WIDTH/MINIMUM_CENTER
    units = 0
    for i in range(1, PANELS+1):
        y = CORE*i/PANELS
        shape = max(F(0), exp_lower(-4*y**2, 25))
        action = DEPTH*minimum_scaled_t**2*shape/2 - (y+HALF_WIDTH)**2/2
        entry = exp_lower(action, 11)
        units += entry.numerator*SCALE//entry.denominator
    ratio = NORMAL_LOWER*2*CORE/PANELS*F(units,SCALE)
    assert ratio == F(214918198851,214748364800) > 1
    print(json.dumps({'aAtLeast':str(MINIMUM_CENTER),'depth':str(DEPTH),'width':'1/2','power':2,
        'coreHalfWidth':str(CORE),'observationHalfWidth':str(HALF_WIDTH),'panelsPerHalf':PANELS,
        'strictIntegratedRatioLower':str(ratio)},indent=2))


if __name__=='__main__': main()
