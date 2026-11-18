"""Analytic five-switch half-circle and certified all-pose obstruction on a continuous tangent band."""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import platform
import sys
import flint
from flint import arb, arb_series, ctx

parser = argparse.ArgumentParser()
parser.add_argument("computation", type=Path)
parser.add_argument("output", type=Path)
args = parser.parse_args()
comp = args.computation.resolve()
src = comp / "b1-tangent-curvature"
sys.path.insert(0, str(src))
from certify_b1_classification import (
    Geometry,
    PINS,
    SELECT,
    dot,
    perp,
    rotate,
    unit,
    rat,
)

expected_base = "584e55bf2b54eaf3273f1d77c03640a8294ae3af9b1b3e1a08daf75f36cbb4dd"
assert (
    sha256((src / "certify_b1_classification.py").read_bytes()).hexdigest()
    == expected_base
)
G = Geometry(comp)
ctx.prec = 192
ctx.cap = 3
pi = arb.pi()
L = pi / 3
p = rat(Q("0.1391310926383"))
X = rat(Q("626.5943632307"))
x2 = X * pi / 2880
x1 = L - ((L - x2).cos() - arb(1) / 2 + L.cos() / 2 - p * L.sin()).acos()
x5 = 2 * L + (arb(3) / 4 + p * L.sin()).acos()


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def scale(a, r):
    return tuple(x * r for x in a)


D = sub(sub(unit(x2), unit(x1)), unit(x5))
m = arb.atan2(D[0], -D[1])
ell = (dot(D, D).sqrt() / 2).asin()
x3 = m - ell
x4 = m + ell
xs = [arb(0), x1, x2, x3, x4, x5, pi]
assert all(a < b for a, b in zip(xs, xs[1:]))
assert x2 < L < x3 < x4 < 2 * L < x5
start = (arb(1) / 2, p)
pieces = []
point = start
for i, (a, b) in enumerate(zip(xs, xs[1:])):
    r = arb(int(i % 2 == 0))
    center = sub(point, scale(unit(a), r))
    pieces.append((r, center))
    point = add(center, scale(unit(b), r))
assert all(v.contains(0) for v in sub(point, sub(start, unit(arb(0)))))


def idx(t):
    indices = [i for i, (a, b) in enumerate(zip(xs, xs[1:])) if a < t < b]
    assert len(indices) == 1, (t, indices)
    return indices[0]


def halfh(t):
    r, center = pieces[idx(t)]
    return r + dot(center, unit(t))


def support(u):
    t = arb.atan2(u[1], u[0])
    assert t > 0 or t < 0
    return halfh(t) if t > 0 else 1 - halfh(t + pi)


assert (halfh(L) - arb(1) / 2).contains(0) and (halfh(2 * L) - arb(1) / 2).contains(0)

# Exact critical-point maximization on each smooth trigonometric piece of F.
events = sorted(
    [arb(0), x1, x2, x3 - L, x4 - L, x5 - 2 * L, L], key=lambda z: float(z.mid())
)
assert len(events) == 7 and all(a < b for a, b in zip(events, events[1:]))
rows = []
for i, (a, b) in enumerate(zip(events, events[1:])):
    C = -arb(1) / 2
    vec = (arb(0), arb(0))
    for j in range(3):
        r, c = pieces[idx((a + b) / 2 + j * L)]
        sign = (-1) ** j
        C += sign * r
        vec = add(vec, scale(rotate(c, -j * L), sign))
    f = lambda t: C + dot(vec, unit(t))
    df = lambda t: dot(vec, perp(unit(t)))
    da, db = df(a), df(b)
    assert (da > 0 or da < 0) and (db > 0 or db < 0)
    ends = []
    if i > 0:
        ends.append(f(a))
    else:
        assert f(a).contains(0) and da < 0
    if i + 1 < len(events) - 1:
        ends.append(f(b))
    else:
        assert f(b).contains(0) and db > 0
    interior = da > 0 and db < 0
    if interior:
        ends.append(C + dot(vec, vec).sqrt())
    upper = max(z.upper() for z in ends)
    assert upper < 0
    rows.append(
        {
            "ends_in_degrees": [str(a * 180 / pi), str(b * 180 / pi)],
            "interior_maximum": interior,
            "nonzero_candidates_upper": str(upper),
            "left_derivative": str(da),
            "right_derivative": str(db),
        }
    )

ah = support(G.normal) - G.height
assert ah > 0
sigma = rat(Q(15497, 10000)) * pi / 180
caps = []
for k in range(0, 360, 60):
    for sign in (1, -1):
        if (k, sign) in ((0, 1), (240, -1)):
            continue
        choice = 3 if SELECT[len(caps)] == 0 else 5
        v = unit(pi / 6 + sigma + choice * pi / 3)
        u = rotate((v[0], sign * v[1]), -k * pi / 180)
        gap = support(u) - arb(1) / 2
        assert gap > 0
        caps.append(
            {"k": k, "reflection_sign": sign, "choice": choice, "margin": str(gap)}
        )

# The vertex after the first material arc is one and the same point for every w.
vertex = G.qa(add(sub(start, unit(arb(0))), unit(x1)))


def jets(w):
    a, b = G.ns["points"](arb_series([w, 1, 0]))[0][2]
    q = G.qa((a, b))
    pos = q[0][0], q[1][0]
    v = q[0][1], q[1][1]
    acc = 2 * q[0][2], 2 * q[1][2]
    speed = dot(v, v).sqrt()
    assert speed > 0
    u = scale(perp(v), -1 / speed)
    du = scale(perp(sub(scale(acc, 1 / speed), scale(v, dot(v, acc) / speed**3))), -1)
    return pos, u, du


lo, hi = Q("0.4157"), Q("0.4420")
N = 512
work = [(lo + (hi - lo) * i / N, lo + (hi - lo) * (i + 1) / N, 0) for i in range(N)]
worst = None
nodes = 0
maxdepth = 0
while work:
    a, b, d = work.pop()
    nodes += 1
    maxdepth = max(maxdepth, d)
    mid = (a + b) / 2
    pos, u, _ = jets(rat(mid))
    gap = dot(sub(vertex, pos), u)
    box, _, du = jets(rat(a).union(rat(b)))
    rate = abs(dot(sub(vertex, box), du)).upper()
    bound = (gap - rat((b - a) / 2) * rate).lower()
    if bound > 0:
        worst = bound if worst is None or bound < worst else worst
    else:
        assert d < 18, (a, b, bound)
        work.extend([(a, mid, d + 1), (mid, b, d + 1)])

out = {
    "python": platform.python_version(),
    "python_flint": flint.__version__,
    "precision_bits": ctx.prec,
    "script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
    "upstream_base_sha256": sha256(
        (src / "certify_b1_classification.py").read_bytes()
    ).hexdigest(),
    "geometry_inputs_sha256": PINS,
    "parameter_p": "0.1391310926383",
    "x2_in_pi_over_2880": "626.5943632307",
    "switches_in_pi_over_2880": [str(t * 2880 / pi) for t in xs[1:-1]],
    "F_pieces": rows,
    "AH_margin": str(ah),
    "ten_other_pose_margins": caps,
    "tangent_band": {
        "lo": str(lo),
        "hi": str(hi),
        "nodes": nodes,
        "leaves": (nodes + N) // 2,
        "initial_cells": N,
        "max_depth": maxdepth,
        "lower_margin": str(worst),
        "fixed_QA_vertex": [str(y) for y in vertex],
    },
}
args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
print(
    "five arcs, six F pieces, AH",
    ah,
    "tangents",
    lo,
    hi,
    "nodes",
    nodes,
    "margin",
    worst,
)
