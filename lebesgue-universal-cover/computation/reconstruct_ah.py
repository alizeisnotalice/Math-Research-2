#!/usr/bin/env python3
"""Numerical coordinate reconstruction of Gibbs (2018), Definition 8.

This is deliberately a geometry audit, not a verification of Proposition 7.
It fixes the directed-angle convention visible in Figure 7: ED is directed
from E1 to D1 and ``s`` is positive clockwise from this direction.  With that
convention L(-theta)=G exactly (up to roundoff).  Gibbs's text does not state
this orientation algebraically, so the convention is recorded rather than
silently assumed.
"""
from __future__ import annotations

from math import atan2, cos, pi, radians, sin, sqrt

from reconstruct_s_primitives import (
    TOL, add, cross, dot, mul, norm, primitives, sub,
    circle_circle_intersections,
)


def rot(a, angle):
    c, s = cos(angle), sin(angle)
    return (c * a[0] - s * a[1], s * a[0] + c * a[1])


def solve_2x2_rows(a, b, ra, rb):
    """Solve a.x=ra and b.x=rb."""
    det = cross(a, b)
    assert abs(det) > TOL
    return ((ra * b[1] - a[1] * rb) / det,
            (a[0] * rb - ra * b[0]) / det)


def angle_from_to(a, b):
    return atan2(cross(a, b), dot(a, b))


def point_in_polygon(x, vertices):
    """Odd-even test for a *sampled* boundary; numerical diagnostic only."""
    inside = False
    xx, yy = x
    for a, b in zip(vertices, vertices[1:] + vertices[:1]):
        ay, by = a[1], b[1]
        if (ay > yy) != (by > yy):
            xcross = a[0] + (yy - ay) * (b[0] - a[0]) / (by - ay)
            if xx < xcross:
                inside = not inside
    return inside


def definition8(sigma_deg):
    """Return a concrete realization of the Figure-7/Definition-8 data."""
    outer, d, q = primitives(sigma_deg)

    # Outer side normals in reconstruct_s_primitives: n_5 is DC, n_3 is EF.
    n3 = (-1.0, 0.0)
    n5 = (0.5, -sqrt(3.0) / 2.0)
    E1, D1 = outer["E1"], outer["D1"]
    M = mul(0.5, add(E1, D1))

    # At s=-theta the unit segment is F3--G.  Figure 7 puts theta between
    # that directed segment and ED.  ED is E1 -> D1; standard angles are CCW,
    # whereas the paper's signed s is taken clockwise in this realization.
    e_ed = mul(1.0 / norm(sub(D1, E1)), sub(D1, E1))
    d_g = mul(1.0 / norm(sub(q["G"], d["F3"])), sub(q["G"], d["F3"]))
    theta = angle_from_to(e_ed, d_g)
    assert theta > 0.0

    # N(t) is on DE, and the supplied endpoint condition N(tau)=E3 determines
    # both the sign and tau; N(-tau) need not be a named dodecagon vertex.
    e_n = mul(1.0 / norm(sub(d["E3"], M)), sub(d["E3"], M))
    tau = norm(sub(d["E3"], M))

    def L(s):
        v = rot(e_ed, -s)                 # unit vector U(s) -> L(s)
        # L lies on DC: n5.L=1/2; L-v lies on EF: n3.(L-v)=1/2.
        return solve_2x2_rows(n5, n3, 0.5, 0.5 + dot(n3, v))

    def U(s):
        v = rot(e_ed, -s)
        return sub(L(s), v)

    def N(t):
        return add(M, mul(t, e_n))

    def X(s, t):
        # Definition 8's "near the top corner A" selects the upper circle
        # crossing.  Its uniqueness must be checked over the parameter box.
        roots = circle_circle_intersections(L(s), N(t))
        return max(roots, key=lambda x: x[1])

    P = X(-theta, -tau)
    X_as = X(-theta, tau)                 # the X already used in Definition 6
    Y = X(theta, tau)                     # source of the shortest path
    Z = X(0.0, 0.0)                       # destination of the shortest path
    W = X(theta, -tau)

    # Boundary parametrization of R, following the rectangle boundary.
    k = 240
    left = [X(-theta, -tau + 2 * tau * j / k) for j in range(k + 1)]
    top = [X(-theta + 2 * theta * j / k, tau) for j in range(k + 1)]
    right = [X(theta, tau - 2 * tau * j / k) for j in range(k + 1)]
    bottom = [X(theta - 2 * theta * j / k, -tau) for j in range(k + 1)]
    boundary = left + top[1:] + right[1:] + bottom[1:]

    # This samples the direct candidate [Y,Z] against that boundary.  It is
    # explicitly not a proof that the straight segment lies in R.
    direct_samples = [add(Y, mul(j / 1000, sub(Z, Y))) for j in range(1001)]
    direct_sample_inside = sum(point_in_polygon(x, boundary) for x in direct_samples)

    return {
        "outer": outer, "d": d, "q": q, "M": M, "theta": theta,
        "tau": tau, "e_ed": e_ed, "e_n": e_n, "L": L, "U": U,
        "N": N, "X": X, "P": P, "X_as": X_as, "Y": Y, "Z": Z,
        "W": W, "boundary": boundary, "direct_sample_inside": direct_sample_inside,
    }


def show(name, x):
    print(f"{name:>20s} = ({x[0]: .15f}, {x[1]: .15f})")


if __name__ == "__main__":
    sigma = 1.5494
    a = definition8(sigma)
    print(f"sigma_deg={sigma:.4f}")
    print(f"theta_deg={a['theta'] * 180/pi:.15f}")
    print(f"tau={a['tau']:.15f}")
    for name in ("M", "e_ed", "e_n", "P", "X_as", "Y", "Z", "W"):
        show(name, a[name])
    for s, t, label in ((-a["theta"], -a["tau"], "P"),
                        (-a["theta"], a["tau"], "X_as"),
                        (a["theta"], a["tau"], "Y"),
                        (a["theta"], -a["tau"], "W"),
                        (0.0, 0.0, "Z")):
        ell, nn, xx = a["L"](s), a["N"](t), a["X"](s, t)
        assert abs(norm(sub(xx, ell)) - 1.0) < 3e-13
        assert abs(norm(sub(xx, nn)) - 1.0) < 3e-13
        print(f"  X({s*180/pi:+.9f} deg, {t:+.12f})={label}; "
              f"|X-L|={norm(sub(xx,ell)):.15f}; |X-N|={norm(sub(xx,nn)):.15f}")
    assert norm(sub(a["L"](-a["theta"]), a["q"]["G"])) < 3e-13
    assert norm(sub(a["U"](-a["theta"]), a["d"]["F3"])) < 3e-13
    assert norm(sub(a["N"](a["tau"]), a["d"]["E3"])) < 3e-13
    # The source defines the negative endpoint by signed distance from M,
    # not as a named dodecagon vertex.  At nonzero slant it is not D2.
    assert norm(sub(a["P"], a["q"]["P"])) < 3e-13
    assert norm(sub(a["X_as"], a["q"]["X"])) < 3e-13
    distances = []
    for i in range(201):
        s = -a["theta"] + 2 * a["theta"] * i / 200
        for j in range(201):
            t = -a["tau"] + 2 * a["tau"] * j / 200
            distances.append(norm(sub(a["N"](t), a["L"](s))))
    print(f"sampled |N(t)-L(s)| range on 201x201 box: "
          f"[{min(distances):.15f}, {max(distances):.15f}]")
    print(f"|Y-G|={norm(sub(a['Y'], a['q']['G'])):.15f}")
    print("direct Y->Z sampled-in-R count "
          f"{a['direct_sample_inside']}/1001 (polygonal diagnostic only)")
