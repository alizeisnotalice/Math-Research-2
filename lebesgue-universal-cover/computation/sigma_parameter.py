"""A shared exact degree interval for the parameter-family replay."""
from fractions import Fraction
import os

from rational_trig import PI, mul

spec = os.environ.get("LEBESGUE_SIGMA_DEGREES", "1.5494,1.5494")
SIGMA = tuple(Fraction(value) for value in spec.split(","))
assert len(SIGMA) == 2 and 0 < SIGMA[0] <= SIGMA[1] < 30
SIGMA_OUT = [str(value) for value in SIGMA]


def radians(base, coefficient=0):
    """Enclose pi*base + coefficient*sigma, with input sigma in degrees."""
    endpoints = [base + coefficient * value / 180 for value in SIGMA]
    return mul(PI, (min(endpoints), max(endpoints)))
