#!/usr/bin/env python3
"""Reconstruct Definition 3 of Gibbs (2018) in fixed coordinates.

The regular hexagon has inradius 1/2, with vertices A,...,F in clockwise
order starting at the top.  The two extra supporting half-planes inferred
from Figure 3 cut C and E.  This is only the P(sigma) base cover: it does
not implement the later arc cuts AH and EH.
"""
from __future__ import annotations

from math import cos, pi, radians, sin, sqrt


EPS = 1.0e-12


def cross(a, b):
    return a[0] * b[1] - a[1] * b[0]


def clip(poly, normal):
    """Intersect a counterclockwise polygon with normal.dot(x) <= 1/2."""
    out = []
    def val(p):
        return normal[0] * p[0] + normal[1] * p[1] - 0.5
    for p, q in zip(poly, poly[1:] + poly[:1]):
        vp, vq = val(p), val(q)
        pin, qin = vp <= EPS, vq <= EPS
        if pin:
            out.append(p)
        if pin != qin:
            t = vp / (vp - vq)
            out.append((p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])))
    return out


def signed_area(poly):
    return 0.5 * sum(cross(p, q) for p, q in zip(poly, poly[1:] + poly[:1]))


def p_sigma(sigma_deg):
    s = radians(sigma_deg)
    # A,B,C,D,E,F, counterclockwise for shoelace orientation.
    h = [(cos(pi / 2 + k * pi / 3) / 3**0.5,
          sin(pi / 2 + k * pi / 3) / 3**0.5) for k in range(6)]
    # The list above is CCW: A,F,E,D,C,B; clipping is orientation independent.
    cut_c = (cos(-pi / 6 + s), sin(-pi / 6 + s))
    cut_e = (cos(-5 * pi / 6 + s), sin(-5 * pi / 6 + s))
    return clip(clip(h, cut_c), cut_e), cut_c, cut_e


def p_sigma_area_formula(sigma_deg):
    """Closed form for the Definition-3 polygon, before any arc cuts.

    The C-corner cut has side offsets d and ell on BC and CD.  Rotation by
    -120 degrees maps it to the E-corner cut, so the total removed area is
    twice (sqrt(3)/4)*d*ell.
    """
    beta = -pi / 6 + radians(sigma_deg)
    ux, uy = cos(beta), sin(beta)
    c = (0.5, -1 / (2 * sqrt(3)))
    cd_direction = (-sqrt(3) / 2, -0.5)  # unit vector C -> D
    d = (0.5 - 0.5 * ux) / uy - c[1]
    ell = (0.5 - (ux * c[0] + uy * c[1])) / (
        ux * cd_direction[0] + uy * cd_direction[1]
    )
    return sqrt(3) / 2 - sqrt(3) / 2 * d * ell


if __name__ == "__main__":
    for angle in (0.0, 1.0, 1.5494):
        poly, cut_c, cut_e = p_sigma(angle)
        print(f"sigma_deg={angle:.4f}")
        print(" cut_C_normal=", cut_c, "cut_E_normal=", cut_e)
        clip_area = abs(signed_area(poly))
        formula_area = p_sigma_area_formula(angle)
        print(" vertices=", len(poly), "area=", f"{clip_area:.15f}")
        print(" formula_area=", f"{formula_area:.15f}", "residual=", f"{clip_area-formula_area:.3e}")
        assert abs(clip_area - formula_area) < 1.0e-12
        for x, y in poly:
            print(f"  ({x:.15f}, {y:.15f})")
    assert abs(p_sigma_area_formula(0.0) - (2 - 2 / sqrt(3))) < 1.0e-12
