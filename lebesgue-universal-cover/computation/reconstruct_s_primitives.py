#!/usr/bin/env python3
"""Coordinate reconstruction of Gibbs (2018), Definitions 4--6.

The labels follow the dodecagon convention recorded in Baez--Bagdasaryan--
Gibbs (2015), p. 7: H cap H' has cyclic vertices
A3,A2,B3,B2,C3,C2,D3,D2,E3,E2,F3,F2.  All calculations here are geometric
definitions only; they do not establish the universal-cover propositions.
"""
from __future__ import annotations

from math import atan2, cos, pi, radians, sin, sqrt


TOL = 1.0e-10


def add(a, b): return (a[0] + b[0], a[1] + b[1])
def sub(a, b): return (a[0] - b[0], a[1] - b[1])
def mul(t, a): return (t * a[0], t * a[1])
def dot(a, b): return a[0] * b[0] + a[1] * b[1]
def norm(a): return sqrt(dot(a, a))
def unit(theta): return (cos(theta), sin(theta))
def tangent(n): return (-n[1], n[0])


def intersection_of_support_lines(a, b):
    """Solve a.x=b.x=1/2 for distinct unit normals a,b."""
    det = a[0] * b[1] - a[1] * b[0]
    assert abs(det) > TOL
    return ((b[1] - a[1]) / (2 * det), (a[0] - b[0]) / (2 * det))


def project_to_support_line(c, n):
    """Orthogonal projection onto n.x=1/2."""
    return sub(c, mul(dot(n, c) - 0.5, n))


def circle_line_intersections(c, n):
    """Intersections of |x-c|=1 and n.x=1/2, in the tangent basis."""
    t = tangent(n)
    q0 = dot(c, t)
    d = dot(c, n) - 0.5
    h2 = 1 - d * d
    assert h2 >= -TOL, h2
    h = sqrt(max(0, h2))
    base = mul(0.5, n)
    return [add(base, mul(q0 + sign * h, t)) for sign in (-1, 1)]


def circle_circle_intersections(a, b):
    """Intersections of the two unit circles centered at a,b."""
    dvec = sub(b, a)
    d = norm(dvec)
    assert 0 < d <= 2 + TOL, d
    mid = mul(0.5, add(a, b))
    perp = (-dvec[1] / d, dvec[0] / d)
    h = sqrt(max(0, 1 - d * d / 4))
    return [add(mid, mul(sign * h, perp)) for sign in (-1, 1)]


def cross(a, b): return a[0] * b[1] - a[1] * b[0]


def arc_green_integral(a, b, center):
    """Integral x dy-y dx along the minor unit-circle arc a -> b."""
    va, vb = sub(a, center), sub(b, center)
    angle = atan2(cross(va, vb), dot(va, vb))
    return cross(center, sub(b, a)) + angle


def curvilinear_area(segments, arcs):
    """Absolute Green-formula area for an oriented piecewise curve."""
    integral = sum(cross(a, b) for a, b in segments)
    integral += sum(arc_green_integral(a, b, c) for a, b, c in arcs)
    return abs(0.5 * integral)


def on_segment(x, a, b):
    ab = sub(b, a)
    q = dot(sub(x, a), ab) / dot(ab, ab)
    return -TOL <= q <= 1 + TOL, q


def labels(sigma_deg):
    s = radians(sigma_deg)
    n = [unit(j * pi / 3) for j in range(6)]
    m = [unit(pi / 6 + s + j * pi / 3) for j in range(6)]
    # Outer regular-hexagon vertices, as in Figure 1 of BBG 2015.
    outer = {
        "A1": intersection_of_support_lines(n[1], n[2]),
        "B1": intersection_of_support_lines(n[0], n[1]),
        "C1": intersection_of_support_lines(n[0], n[5]),
        "D1": intersection_of_support_lines(n[4], n[5]),
        "E1": intersection_of_support_lines(n[3], n[4]),
        "F1": intersection_of_support_lines(n[2], n[3]),
    }
    # The 12 vertices of H cap H', in the BBG/Gibbs naming convention.
    d = {
        "A3": intersection_of_support_lines(n[2], m[1]),
        "A2": intersection_of_support_lines(n[1], m[1]),
        "B3": intersection_of_support_lines(n[1], m[0]),
        "B2": intersection_of_support_lines(n[0], m[0]),
        "C3": intersection_of_support_lines(n[0], m[5]),
        "C2": intersection_of_support_lines(n[5], m[5]),
        "D3": intersection_of_support_lines(n[5], m[4]),
        "D2": intersection_of_support_lines(n[4], m[4]),
        "E3": intersection_of_support_lines(n[4], m[3]),
        "E2": intersection_of_support_lines(n[3], m[3]),
        "F3": intersection_of_support_lines(n[3], m[2]),
        "F2": intersection_of_support_lines(n[2], m[2]),
    }
    return n, m, outer, d


def primitives(sigma_deg):
    n, m, outer, d = labels(sigma_deg)
    # Definition 4: K is tangency on C3C2 (m[5]); G is the root on D3C2 (n[5]).
    K = project_to_support_line(d["F3"], m[5])
    G_roots = circle_line_intersections(d["F3"], n[5])
    G_candidates = [(x, on_segment(x, d["D3"], d["C2"])[1]) for x in G_roots
                    if on_segment(x, d["D3"], d["C2"])[0]]
    assert len(G_candidates) == 1, G_candidates
    G = G_candidates[0][0]

    # Definition 5: tangencies I,H then the two possible unit-circle crossings J.
    I = project_to_support_line(d["B3"], m[3])
    H = project_to_support_line(d["C3"], n[3])
    J_roots = circle_circle_intersections(d["B3"], d["C3"])

    # Definition 6: P,Q tangencies and the two candidate circle crossings X.
    P = project_to_support_line(G, n[2])
    Q = project_to_support_line(d["E3"], n[1])
    X_roots = circle_circle_intersections(G, d["E3"])
    # In the pictured branch the desired crossings are respectively nearest
    # E2 and A1.  The other root is visibly outside the base hexagon.
    J = min(J_roots, key=lambda x: norm(sub(x, d["E2"])))
    X = min(X_roots, key=lambda x: norm(sub(x, outer["A1"])))

    areas = {
        # K -> G follows the F3 arc; the remaining two pieces are segments.
        "CS": curvilinear_area([(G, d["C2"]), (d["C2"], K)],
                                [(K, G, d["F3"])]),
        # I -> J is centered at B3; J -> H is centered at C3.
        "ES": curvilinear_area([(H, d["E2"]), (d["E2"], I)],
                                [(I, J, d["B3"]), (J, H, d["C3"])]),
        # Q -> X is centered at E3; X -> P is centered at G.
        "AS": curvilinear_area([(P, outer["A1"]), (outer["A1"], Q)],
                                [(Q, X, d["E3"]), (X, P, G)]),
    }
    return outer, d, {"K": K, "G": G, "I": I, "H": H,
                       "J": J, "J_roots": J_roots, "P": P, "Q": Q,
                       "X": X, "X_roots": X_roots, "areas": areas}


def show(name, x):
    print(f"{name:>8s} = ({x[0]: .15f}, {x[1]: .15f})")


if __name__ == "__main__":
    for sigma in (0.0, 1.3, 1.5494):
        print(f"sigma_deg={sigma:.4f}")
        outer, d, p = primitives(sigma)
        for name in ("F3", "C3", "C2", "D3", "B3", "E3", "E2"):
            show(name, d[name])
        for name in ("K", "G", "I", "H", "J", "P", "Q", "X"):
            show(name, p[name])
        for name in ("J_roots", "X_roots"):
            for k, x in enumerate(p[name]): show(f"{name}[{k}]", x)
        # Defining equalities and the intended line-segment branch for G.
        assert abs(norm(sub(p["K"], d["F3"])) - 1) < TOL
        assert abs(norm(sub(p["G"], d["F3"])) - 1) < TOL
        assert abs(norm(sub(p["I"], d["B3"])) - 1) < TOL
        assert abs(norm(sub(p["H"], d["C3"])) - 1) < TOL
        assert abs(norm(sub(p["P"], p["G"])) - 1) < TOL
        assert abs(norm(sub(p["Q"], d["E3"])) - 1) < TOL
        assert on_segment(p["K"], d["C3"], d["C2"])[0]
        assert on_segment(p["G"], d["D3"], d["C2"])[0]
        assert on_segment(p["I"], d["E3"], d["E2"])[0]
        assert on_segment(p["H"], d["E2"], d["F3"])[0]
        assert on_segment(p["P"], outer["F1"], outer["A1"])[0]
        assert on_segment(p["Q"], outer["A1"], outer["B1"])[0]
        print(" areas=", ", ".join(f"{name}:{value:.15e}" for name, value in p["areas"].items()))
