#!/usr/bin/env python3
"""Exact rational arithmetic for the fixed profile's bounds and subsequential rates."""
from fractions import Fraction as F


def cost(u):
    return u*u/32 + (1-u)**2/2


def rate(c):
    minimizer = F(16, 17)
    u = 1/c
    while u > minimizer:
        u /= 4
    while 4*u <= minimizer:
        u *= 4
    return max(F(-1, 8), F(1, 16)-min(cost(u), cost(4*u)))


assert rate(F(1)) == F(1, 32)
assert rate(F(5, 4)) == F(9, 400)
assert rate(F(2)) == F(-9, 128)

plateau_log_lower = F(320**2, 16) - F(256**2, 32) - (64+F(1, 16))**2/2 - 4
assert plateau_log_lower > 2295
assert F(64)-F(1, 8) > 1

left_log_upper = F(512**2, 16)-(256-F(1, 8))**2/2
right_log_upper = F(512**2, 16)-F(1024**2, 32)
assert max(left_log_upper, right_log_upper, F(-512**2, 8)) < -15000
print("Exact full rates at c=1,5/4,2", *(str(rate(c)) for c in (F(1),F(5,4),F(2))))
print("Exact unit-window rate at c=5/4", F(-1,8))
print("λ=320: log I(full) >", plateau_log_lower, "> 2295; I(unit)<1")
print("λ=512: I(full) < 6 exp(−15000) < 1")
