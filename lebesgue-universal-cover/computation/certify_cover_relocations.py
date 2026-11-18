"""Check exact normal maps, exhaustive case logic and fixed-angle line orders."""
import contextlib
from fractions import Fraction as F
import io
import itertools
import json
from pathlib import Path
import runpy
import sys

with contextlib.redirect_stdout(io.StringIO()):
    base = runpy.run_path(str(Path(__file__).with_name("certify_ah_path.py")))
I, V, add, sub, mul, dot, perp = (base[name] for name in ("I", "V", "add", "sub", "mul", "dot", "perp"))
PI, scale, trig = (base[name] for name in ("PI", "scale", "trig"))
half = I(F(1, 2))
sigma = F("1.5494") / 180


def unit(q):
    angle = scale(PI, q)
    return V(I(*trig(angle, "cos")), I(*trig(angle, "sin")))


def cross(a, b):
    return a[0] * b[1] - a[1] * b[0]


def meet(a, b, ra=half, rb=half):
    determinant = cross(a, b)
    return V((ra * b[1] - a[1] * rb) / determinant,
             (a[0] * rb - ra * b[0]) / determinant)


def signed(value, direction):
    assert (value.lo > 0 if direction == 1 else value.hi < 0), value.out()
    return value.out()


n = [unit(F(j, 3)) for j in range(6)]
pos = [unit(F(1, 6) + sigma + F(j, 3)) for j in range(6)]
neg = [unit(F(1, 6) - sigma + F(j, 3)) for j in range(6)]

# The F1/C1 reflection exchanges the two F cuts and fixes their intersection.
R = meet(pos[2], neg[2])
F3p, F2p = meet(n[3], pos[2]), meet(n[2], pos[2])
F3m, F2m = meet(n[3], neg[2]), meet(n[2], neg[2])
C3 = meet(n[0], pos[5])
# Exact S_x=1/2 and S_y below R_y determine the Definition-9 root.
S = V(half, R[1] - (1 - (half - R[0]) * (half - R[0])).sqrt())

checks = {
    "R_is_strictly_inside_original_hexagon": [
        signed(half - dot(normal, R), 1) for normal in n],
    "C3_below_S": signed(S[1] - C3[1], 1),
    "R_above_S": signed(R[1] - S[1], 1),
    "C3_voronoi": signed(dot(sub(C3, R), sub(F3p, R)), -1),
    "S_voronoi": signed(dot(sub(S, R), sub(F3p, R)), -1),
    "negative_cut_at_F3plus": signed(dot(neg[2], F3p) - half, -1),
    "negative_cut_at_F2plus": signed(dot(neg[2], F2p) - half, 1),
    "positive_cut_at_F3minus": signed(dot(pos[2], F3m) - half, 1),
    "positive_cut_at_F2minus": signed(dot(pos[2], F2m) - half, -1),
    "positive_F_normal_y": signed(pos[2][1], 1),
    "negative_F_normal_y": signed(neg[2][1], 1),
    "F3plus_above_S": signed(F3p[1] - S[1], 1),
    "F3minus_above_S": signed(F3m[1] - S[1], 1),
    "S_above_EH_slab": signed(S[1] + F(133,1000), 1),
}

# Angles have units pi/6, so m_j^e is encoded by (1+2j,e) mod 12,
# with coefficient e of sigma. Original normals are (2j,0).
def angle(j, e):
    return ((1 + 2 * j) % 12, e)


def reflect(axis, direction):
    k, e = direction
    return ((2 * axis - k) % 12, -e)


def rotate(k, direction):
    a, e = direction
    return ((a + k) % 12, e)


ah_map = [(4, 3), (2, 5)]
q1_map = [(5, 5), (3, 1)]
q2_map = [(5, 4), (3, 0)]
for axis, mapping in [(2, ah_map), (5, q1_map), (4, q2_map)]:
    for plus, minus in mapping:
        if axis == 2:
            assert reflect(axis, angle(plus, -1)) == angle(minus, 1)
        else:
            assert reflect(axis, angle(plus, 1)) == angle(minus, -1)
assert reflect(5, (4, 0)) == (6, 0)  # Q1 maps F2- (n2) onto F3+ (n3).
assert reflect(5, angle(2, -1)) == angle(2, 1)
assert reflect(4, (1, 0)) == (7, 0)  # Q2 reverses the EH strip normal.
for j, target, sign in [(3, 5, 1), (1, 3, 1), (2, 4, -1), (0, 2, -1)]:
    assert rotate(4, angle(j, sign)) == angle(target, sign)

count = 0
for A, D, Ap, Dp in itertools.product((False, True), repeat=4):
    if (A and D) or (Ap and Dp):
        continue
    assert (not Ap) or (not D and not Dp) or (not A and not Dp)
    count += 1

print(json.dumps({"python": sys.version, "sigma_degrees": "7747/5000",
                  "consistent_hit_patterns": count, "exact_maps_passed": True,
                  "fixed_angle_intervals": checks}, indent=2))
