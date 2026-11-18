#!/usr/bin/env python3
"""Rational-interval audit of the fixed-angle Gibbs Definition-9 EH splice.

This is a fixed-sigma geometric certificate.  It proves neither Gibbs's
orbiform relocation argument nor universal coverage.
"""
from __future__ import annotations

import contextlib
from fractions import Fraction as F
import hashlib
import io
import json
from pathlib import Path
import runpy
import sys

base = Path(__file__).parent
with contextlib.redirect_stdout(io.StringIO()):
    d = runpy.run_path(str(base / "certify_ah_path.py"))

I = d["I"]
PI = d["PI"]
trig = d["trig"]
scale = d["scale"]
add, sub, mul, dot, norm2, perp = (d[k] for k in ("add", "sub", "mul", "dot", "norm2", "perp"))
r3 = d["r3"]


def cross(a, b):
    return a[0] * b[1] - a[1] * b[0]


def unit(q: F):
    angle = scale(PI, q)
    return (I(*trig(angle, "cos")), I(*trig(angle, "sin")))


def meet(a, b, ra=I(F(1, 2)), rb=I(F(1, 2))):
    det = cross(a, b)
    return ((ra * b[1] - a[1] * rb) / det, (a[0] * rb - ra * b[0]) / det)


def project(c, n):
    return sub(c, mul(dot(n, c) - F(1, 2), n))


def circle_line(c, n, sign: int):
    distance = dot(n, c) - F(1, 2)
    return add(project(c, n), mul(sign * (1 - distance * distance).sqrt(), perp(n)))


def circle_circle(a, b, sign: int):
    q = sub(b, a)
    dsq = norm2(q)
    return add(mul(I(F(1, 2)), add(a, b)),
               mul(sign * (1 - dsq / 4).sqrt() / dsq.sqrt(), perp(q)))


def tangent_ccw(point, centre):
    """Unit tangent to the counterclockwise arc at ``point``."""
    return perp(sub(point, centre))


def line_intersection_parameters(a, b, c, e):
    """a+s(b-a)=c+t(e-c), with interval-certified nonzero determinant."""
    q, r = sub(b, a), sub(e, c)
    det = cross(q, r)
    assert det.lo > 0 or det.hi < 0
    return cross(sub(c, a), r) / det, cross(sub(c, a), q) / det


def out(x):
    return x.out()


sigma = F("1.5494")
n = [unit(F(j, 3)) for j in range(6)]
m = [unit(F(1, 6) + sigma / 180 + F(j, 3)) for j in range(6)]
B3 = meet(n[1], m[0])
C3 = meet(n[0], m[5])
E2 = meet(n[3], m[3])

# Definition 5 endpoints, retained to audit which portions of ES are replaced.
Ip = project(B3, m[3])
H = project(C3, n[3])
J = circle_circle(B3, C3, -1)
J_other = circle_circle(B3, C3, 1)

# Definition 9 endpoints, directly from its supporting lines/circles.
R = meet(m[2], n[1], I(F(1, 2)), I(0))
S = circle_line(R, n[0], -1)
S_other = circle_line(R, n[0], 1)
T = circle_circle(S, B3, 1)
T_other = circle_circle(S, B3, -1)
U = circle_line(C3, unit(F(7, 6)), -1)
U_other = circle_line(C3, unit(F(7, 6)), 1)

# Root identities and visually prescribed selections.  Bounds are all strict.
assert J[0].hi < 0 and J_other[0].lo > 1
assert S[1].hi < -F(1, 10) and S_other[1].lo > F(1, 10)
assert T[0].hi < 0 and T_other[0].lo > 1
assert U[1].lo > -F(1, 5) and U_other[1].hi < -F(9, 10)

q = sub(U, T)                         # T -> U
departure = dot(sub(T, S), q)
assert departure.lo > 0
US_distance_squared = norm2(sub(U, S))
assert US_distance_squared.lo > 1

# The chord meets each of the B3 and C3 circles at its named endpoint only.
# Parameter is T + t(U-T); the second B3 root is tB, and the other C3 root tC.
q2 = norm2(q)
tB = -2 * dot(sub(T, B3), q) / q2
tC = -2 * dot(sub(T, C3), q) / q2 - 1
tS = -2 * dot(sub(T, S), q) / q2
assert tB.lo > 1 and tC.hi < 0 and tS.hi < 0

# On the S-circle, |T+t(U-T)-S|^2-1 = 2t*departure+t^2|q|^2.
# Thus the chord is strictly outside the S disk away from T.

# T lies strictly in the short I--J B3 arc; U lies strictly in H--J C3 arc.
def arc_signs(a, b, c):
    return cross(sub(a, c), sub(b, c)), dot(sub(a, c), sub(b, c))

IJ_cross, IJ_dot = arc_signs(Ip, J, B3)
IT_cross, IT_dot = arc_signs(Ip, T, B3)
TJ_cross, TJ_dot = arc_signs(T, J, B3)
HJ_cross, HJ_dot = arc_signs(H, J, C3)
HU_cross, HU_dot = arc_signs(H, U, C3)
UJ_cross, UJ_dot = arc_signs(U, J, C3)
assert IJ_cross.hi < 0 and IT_cross.hi < 0 and TJ_cross.hi < 0
assert IJ_dot.lo > 0 and IT_dot.lo > 0 and TJ_dot.lo > 0
assert HJ_cross.lo > 0 and HU_cross.lo > 0 and UJ_cross.lo > 0
assert HJ_dot.lo > 0 and HU_dot.lo > 0 and UJ_dot.lo > 0
# The alternative circle intersection is on neither short arc.
assert arc_signs(Ip, J_other, B3)[1].hi < 0
assert arc_signs(U, J_other, C3)[1].hi < 0

# The only potentially nearby non-arc parts of ES are H--E2 and E2--I.
# Their line intersections with TU fall outside the two closed segments.
s_HE2, u_HE2 = line_intersection_parameters(T, U, H, E2)
s_E2I, u_E2I = line_intersection_parameters(T, U, E2, Ip)
assert s_HE2.lo > 1 and s_E2I.hi < 0

# Positive orientation of the surviving splice H -> U -> T -> I.
# The two circular portions have positive curvature, and both new corners turn
# strictly left through an angle less than pi/2 (cross and dot are positive).
in_U = tangent_ccw(U, C3)              # H -> U
chord_UT = mul(-1, q)                   # U -> T
out_T = tangent_ccw(T, B3)              # T -> I
turn_U_cross, turn_U_dot = cross(in_U, chord_UT), dot(in_U, chord_UT)
turn_T_cross, turn_T_dot = cross(chord_UT, out_T), dot(chord_UT, out_T)
assert turn_U_cross.lo > 0 and turn_U_dot.lo > 0
assert turn_T_cross.lo > 0 and turn_T_dot.lo > 0
assert HU_cross.lo > 0 and arc_signs(T, Ip, B3)[0].lo > 0

# Definition 9's zero-slant Pál line is the line pal.x=1/2, with
# pal=(cos(7pi/6),sin(7pi/6)).  EH lies on its retained side.  On U--J the
# derivative of pal.x along the CCW arc is cross(x-C3,pal)>0; on J--T it is
# cross(x-B3,pal)<0.  The minor-arc dot tests above rule out a hidden critical
# point between either pair of endpoints.
pal = unit(F(7, 6))
pal_U, pal_J, pal_T = dot(pal, U), dot(pal, J), dot(pal, T)
pal_C_U, pal_C_J = cross(sub(U, C3), pal), cross(sub(J, C3), pal)
pal_B_J, pal_B_T = cross(sub(J, B3), pal), cross(sub(T, B3), pal)
assert pal_U.lo > F(1, 2) - F(1, 10**15)
assert pal_J.lo > F(1, 2) and pal_T.lo > F(1, 2)
assert pal_C_U.lo > 0 and pal_C_J.lo > 0
assert pal_B_J.hi < 0 and pal_B_T.hi < 0

# The same short-arc derivative test shows that both curved EH sides stay
# outside the open unit disk centred at S.  The chord case was proved above.
def outside_derivative(point, centre):
    return 2 * cross(sub(point, centre), sub(centre, S))

S_U = norm2(sub(U, S)) - 1
S_J = norm2(sub(J, S)) - 1
S_T = norm2(sub(T, S)) - 1
S_C_U, S_C_J = outside_derivative(U, C3), outside_derivative(J, C3)
S_B_J, S_B_T = outside_derivative(J, B3), outside_derivative(T, B3)
assert S_U.lo > 0 and S_J.lo > 0 and S_T.lo < 0 < S_T.hi
assert S_C_U.lo > 0 and S_C_J.lo > 0
assert S_B_J.hi < 0 and S_B_T.hi < 0

# E' is the negative-slant E cap: eprime.x>=1/2 for
# eprime=u(7pi/6-sigma).  It contains EH strictly, including all three
# boundary pieces at this fixed positive slant.
eprime = unit(F(7, 6) - sigma / 180)
ep_U, ep_J, ep_T = dot(eprime, U), dot(eprime, J), dot(eprime, T)
ep_C_U, ep_C_J = cross(sub(U, C3), eprime), cross(sub(J, C3), eprime)
ep_B_J, ep_B_T = cross(sub(J, B3), eprime), cross(sub(T, B3), eprime)
assert ep_U.lo > F(1, 2) and ep_J.lo > F(1, 2) and ep_T.lo > F(1, 2)
assert ep_C_U.lo > 0 and ep_C_J.lo > 0
assert ep_B_J.hi < 0 and ep_B_T.hi < 0

# A simple separation audit against the non-neighbouring cuts.  CS is wholly
# in x>0.37 (its K--G F3-arc has negative radial y and clockwise direction),
# while AS and the certified AH boundary are in y>0.49; EH is in x<-0.499 and
# -0.143<y<-0.133.  The numerical constants below are intentionally loose and
# are checked from the exact interval endpoint constructions.
K = project(meet(n[3], m[2]), m[5])
F3 = meet(n[3], m[2])
G = circle_line(F3, n[5], 1)
P = project(G, n[2])
E3 = meet(n[4], m[3])
Q = project(E3, n[1])
X = circle_circle(G, E3, -1)
assert K[0].lo > F(37, 100) and G[0].lo > F(37, 100)
assert P[1].lo > F(49, 100) and Q[1].lo > F(49, 100) and X[1].lo > F(49, 100)
assert U[0].hi < -F(499, 1000) and T[0].hi < -F(499, 1000)
assert T[1].lo > -F(143, 1000) and U[1].hi < -F(133, 1000)
# Monotonicity gives the same bounds on the circular sides, not just endpoints.
CS_arc_cross, _ = arc_signs(K, G, F3)
assert CS_arc_cross.hi < 0 and sub(K, F3)[1].hi < 0 and sub(G, F3)[1].hi < 0
AS_QX_cross, _ = arc_signs(Q, X, E3)
AS_XP_cross, _ = arc_signs(X, P, G)
assert AS_QX_cross.lo > 0 and sub(Q, E3)[0].lo > 0 and sub(X, E3)[0].lo > 0
assert AS_XP_cross.lo > 0 and sub(X, G)[0].hi < 0 and sub(P, G)[0].hi < 0

result = {
    "python": sys.version,
    "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "interval_engine": "rational outward-rounded intervals from certify_ah_path.py",
    "sigma_degrees_exact": str(sigma),
    "points": {k: [out(z) for z in v] for k, v in {
        "B3": B3, "C3": C3, "E2": E2, "I": Ip, "H": H, "J": J,
        "J_other": J_other, "R": R, "S": S, "S_other": S_other,
        "T": T, "T_other": T_other, "U": U, "U_other": U_other,
    }.items()},
    "convex_hull_branch": {
        "departure_dot_T_minus_S__dot__U_minus_T": out(departure),
        "U_distance_to_S_squared": out(US_distance_squared),
        "second_chord_circle_parameters": {"B3": out(tB), "C3": out(tC), "S": out(tS)},
        "conclusion": "departure>0, so TU immediately and thereafter leaves the S disk; the Definition-9 convex-hull boundary is the straight chord TU, not the arc-plus-tangent branch.",
    },
    "arc_order": {
        "I_to_T_to_J_on_B3": {"I_J_cross": out(IJ_cross), "I_T_cross": out(IT_cross), "T_J_cross": out(TJ_cross)},
        "H_to_U_to_J_on_C3": {"H_J_cross": out(HJ_cross), "H_U_cross": out(HU_cross), "U_J_cross": out(UJ_cross)},
        "alternative_intersection_excluded_by_negative_dot": {"B3": out(arc_signs(Ip, J_other, B3)[1]), "C3": out(arc_signs(U, J_other, C3)[1])},
    },
    "ES_compatibility": {
        "shared_arcs": "EH uses U->J and J->T, exactly the reversals of the U->J and J->T subarcs of the positively oriented ES boundary.",
        "TU_vs_HE2_parameters": [out(s_HE2), out(u_HE2)],
        "TU_vs_E2I_parameters": [out(s_E2I), out(u_E2I)],
        "conclusion": "No crossing of the nonshared ES segments; the two common arcs have opposite boundary orientations, hence the two Jordan interiors are locally on opposite sides and have disjoint interiors.",
    },
    "surviving_boundary_positive_turning": {
        "order": "H -> U (C3 arc), U -> T (chord), T -> I (B3 arc)",
        "U": {"cross": out(turn_U_cross), "dot": out(turn_U_dot)},
        "T": {"cross": out(turn_T_cross), "dot": out(turn_T_dot)},
        "conclusion": "both corners are strict left turns of size < pi/2; both retained arcs are counterclockwise unit arcs and hence have positive curvature.",
    },
    "EH_containment": {
        "pal_normal_7pi_over_6_dot": {"U": out(pal_U), "J": out(pal_J), "T": out(pal_T)},
        "pal_arc_derivatives": {"C3_U": out(pal_C_U), "C3_J": out(pal_C_J), "B3_J": out(pal_B_J), "B3_T": out(pal_B_T)},
        "S_disk_squared_distance_minus_one": {"U": out(S_U), "J": out(S_J), "T": out(S_T)},
        "S_disk_arc_derivatives": {"C3_U": out(S_C_U), "C3_J": out(S_C_J), "B3_J": out(S_B_J), "B3_T": out(S_B_T)},
        "Eprime_normal_7pi_over_6_minus_sigma_dot": {"U": out(ep_U), "J": out(ep_J), "T": out(ep_T)},
        "Eprime_arc_derivatives": {"C3_U": out(ep_C_U), "C3_J": out(ep_C_J), "B3_J": out(ep_B_J), "B3_T": out(ep_B_T)},
        "conclusion": "EH is contained in pal.x>=1/2 (equivalently (-pal).x<=-1/2), is strictly contained in the E' cap eprime.x>1/2, and its interior is in the exterior of the S unit disk; equality occurs only at U for the Pál line and at T for the S disk.",
    },
    "non_neighbouring_cut_separation": {
        "CS_x_lower_at_endpoints": [out(K[0]), out(G[0])],
        "AS_y_lower_at_endpoints": [out(P[1]), out(Q[1]), out(X[1])],
        "arc_monotonicity_crosses": {"CS_K_to_G": out(CS_arc_cross), "AS_Q_to_X": out(AS_QX_cross), "AS_X_to_P": out(AS_XP_cross)},
        "EH_endpoint_box": {"x_upper": [out(U[0]), out(T[0])], "y_lower_T": out(T[1]), "y_upper_U": out(U[1])},
        "scope": "Endpoint/radial-monotonicity separation for CS and AS; AH is separately certified in the upper y-region. This is not the global covering proof.",
    },
    "passed": True,
    "scope": "fixed-angle Definition-9 endpoint, convex-hull branch, local Jordan/splice and cut-separation audit only; no universal-cover claim",
}
print(json.dumps(result, indent=2))
