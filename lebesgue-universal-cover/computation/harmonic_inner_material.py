"""Certified common inner material for the fixed seven-body placement of Round 56."""

import json
from flint import arb, ctx
from certify_nine import (
    S,
    PI_QUALITY,
    ar,
    append_point,
    arc_body,
    fourier_body,
    FOURIER_TERMS,
    hull,
)

ctx.prec = PI_QUALITY
DEN = S * S


def polygon(raw, k):
    params = [ar(z) for z in json.loads(raw)["active_parameters"]]
    points = []
    pi = arb.pi()
    for j in range(8 * k):
        s, c = (2 * pi * j / (8 * k)).sin_cos()
        append_point(points, c / 2, s / 2)
    arc_body(points, k, pi, 3, arb(0), *params[:2])
    for n, j in ((5, 2), (7, 5), (9, 8)):
        arc_body(points, k, pi, n, *params[j : j + 3])
    for j, fs in enumerate(FOURIER_TERMS):
        s = 11 + 3 * j
        fourier_body(points, (8 if j == 0 else 16) * k, pi, fs, *params[s : s + 3])
    # Every rounded point is within 2/S of H. Since (1/2)D_unit is in H,
    # scaling the rounded hull by (1-4/S) places it inside H.
    return [(x * (S - 4), y * (S - 4)) for x, y in hull(points)]


def det(a, b):
    return a[0] * b[1] - a[1] * b[0]
