"""Exact variational calibration of the R214 double hit in both remaining tangent gaps."""

import argparse, json, platform, sys
from pathlib import Path
from fractions import Fraction as Q
from hashlib import sha256
import flint
from flint import arb, arb_series, arb_mat, ctx

ap = argparse.ArgumentParser()
ap.add_argument("computation", type=Path)
ap.add_argument("output", type=Path)
args = ap.parse_args()
comp = args.computation.resolve()
src = comp / "b1-tangent-curvature"
sys.path.insert(0, str(src))
from certify_b1_classification import (
    Geometry,
    PINS,
    SELECT,
    rat,
    unit,
    dot,
    perp,
    rotate,
)
from dual_analytic import jets

assert (
    sha256((src / "certify_b1_classification.py").read_bytes()).hexdigest()
    == "584e55bf2b54eaf3273f1d77c03640a8294ae3af9b1b3e1a08daf75f36cbb4dd"
)
assert __debug__
assert (
    sha256((src / "dual_analytic.py").read_bytes()).hexdigest()
    == "1dcc6fed97c14488b84da63b1c76c69662fea01a1362be82f94ddcadaa4aae3b"
)
G = Geometry(comp)
ctx.prec = 224
ctx.cap = 4
pi = arb.pi()
L = pi / 3
step = pi / 2880
fixed = rat(Q("6.1")) * pi / 180
zero = arb(0)
one = arb(1)


def add(u, v):
    return tuple(x + y for x, y in zip(u, v))


def sub(u, v):
    return tuple(x - y for x, y in zip(u, v))


def scale(u, z):
    return tuple(z * x for x in u)


def turn(u, t):
    return (u[0] * t.cos() - u[1] * t.sin(), u[0] * t.sin() + u[1] * t.cos())


def magnitude(u):
    return dot(u, u).sqrt()


def norm_up(u):
    return sum(abs(x).upper() ** 2 for x in u).sqrt().upper()


# The two Jacobian rows use only arithmetic, square roots and trigonometric functions.
def model(p, X):
    x2 = X * step
    e2 = (x2.cos(), x2.sin())
    a = (L - x2).cos() - one / 2 + L.cos() / 2 - p * L.sin()
    b = arb(3) / 4 + p * L.sin()
    e1 = turn((a, -(1 - a * a).sqrt()), L)
    e5 = turn((b, (1 - b * b).sqrt()), 2 * L)
    D = sub(sub(e2, e1), e5)
    r = magnitude(D)
    M = scale(perp(D), (1 - r * r / 4).sqrt() / r)
    e3 = add(M, scale(D, one / 2))
    e4 = sub(M, scale(D, one / 2))
    Ppi = (-one / 2, p)
    A = add(Ppi, e1)
    B = add(sub(A, e2), e3)
    C = sub(Ppi, e5)
    V = add(sub(A, turn(sub(A, e2), -L)), turn(C, -2 * L))
    return [e1, e2, e3, e4, e5], A, V, dot(G.normal, B) - G.height


p_scale = rat(Q("1e-9"))
X_scale = rat(Q("1e-6"))


def critical(kind, p_center, X_center):
    def equations(z):
        _, _, V, d0 = model(p_center + p_scale * z[0], X_center + X_scale * z[1])
        return [
            d0,
            (dot(V, V) - one / 4 if kind == "peak" else dot(V, unit(fixed)) - one / 2),
        ]

    def jac(z):
        columns = []
        for j in range(2):
            columns.append(
                [
                    g[1]
                    for g in equations(
                        [arb_series([u, int(i == j)]) for i, u in enumerate(z)]
                    )
                ]
            )
        return arb_mat([[columns[j][i] for j in range(2)] for i in range(2)])

    box = [arb(0, 1), arb(0, 1)]
    I = arb_mat([[1, 0], [0, 1]])
    contraction = None
    for it in range(4):
        mid = [z.mid() for z in box]
        Y = jac(mid).inv()
        loop = I - Y * jac(box)
        C = arb_mat([[m] for m in mid]) - Y * arb_mat([[z] for z in equations(mid)])
        rem = loop * arb_mat([[box[i] - mid[i]] for i in range(2)])
        new = [C[i, 0] + rem[i, 0] for i in range(2)]
        q = max(sum(abs(loop[i, j]).upper() for j in range(2)) for i in range(2))
        assert q < 1
        if it == 0:
            assert all(box[i].contains_interior(new[i]) for i in range(2))
            contraction = q
        box = new
    p = p_center + p_scale * box[0]
    X = X_center + X_scale * box[1]
    eu, A, V, d0 = model(p, X)
    theta = arb.atan2(V[1], V[0]) if kind == "peak" else fixed
    cuts = [arb.atan2(e[1], e[0]) for e in eu]
    assert zero < cuts[0] < cuts[1] < L < cuts[2] < cuts[3] < 2 * L < cuts[4] < pi
    ah = arb.atan2(G.normal[1], G.normal[0])
    assert cuts[2] < ah < cuts[3]
    assert (
        cuts[0] < theta < cuts[1]
        and cuts[1] < theta + L < cuts[2]
        and cuts[3] < theta + 2 * L < cuts[4]
    )
    return p, X, eu, A, V, theta, cuts, contraction


# Events in a half circle, equality multipliers, and positive parts of the dual kernel.
def sources(theta):
    h = arb.atan2(G.normal[1], G.normal[0])
    assert 0 < h < pi
    return [
        ([(h, 1)], 0, 0),
        ([(theta, -1), (theta + L, 1), (theta + 2 * L, -1)], 0, 0),
        ([(L, 1)], 0, 0),
        ([(2 * L, 1)], 0, 0),
        ([], 1, 0),
        ([], 0, 1),
    ]


def analytic_row(cols, t):
    return [
        sum(
            (c * (a - t).sin() for a, c in events if a > t), d1 * t.sin() + d2 * t.cos()
        )
        for events, d1, d2 in cols
    ]


def calibrate(theta, cuts):
    cols = sources(theta)
    pc = [sum((c * a.sin() for a, c in events), zero) for events, _, _ in cols]
    inv = arb_mat([pc] + [analytic_row(cols, t) for t in cuts]).inv()
    # Q_A^T u's angle is between cut 1 and theta. Its column is sin(tau) in the p-row,
    # and sin(tau-x1) at the first switch, zero at all subsequent switches.
    cosine = inv * arb_mat([[0], [cuts[0].sin()], [0], [0], [0], [0]])
    sine = inv * arb_mat([[-1], [-cuts[0].cos()], [0], [0], [0], [0]])
    # Rigorous checks over all preimage tangent angles in (4.8,5.1) degrees.
    low = rat(Q("4.8")) * pi / 180
    high = rat(Q("5.1")) * pi / 180
    assert cuts[0] < low < high < theta < cuts[1]
    event0 = sources(theta)[0][0][0][0]
    order = [zero, low, theta, L, theta + L, event0, 2 * L, theta + 2 * L, pi]
    assert all(a < b for a, b in zip(order, order[1:]))
    expected = (1, -1, -1, 1, 1, -1, 1, 1, -1)
    mins = None
    for j in range(256):
        tl = low + (high - low) * j / 256
        tr = low + (high - low) * (j + 1) / 256
        tau = tl.union(tr)
        w = [cosine[i, 0] * tau.cos() + sine[i, 0] * tau.sin() for i in range(6)]
        assert w[0] > 0 and w[1] > 0
        events = [zero, tau, theta, L, theta + L, event0, 2 * L, theta + 2 * L, pi]
        vals = []
        assert all(b - a < pi and a < b for a, b in zip(events, events[1:]))
        for e, t in enumerate(events):
            # Skip equal events by their indices. Ordering here is strict on the continuous angle box.
            assert all(x < t or x > t for x in cuts)
            assert expected[e] == (-1) ** sum(x < t for x in cuts)
            amps = analytic_row(cols, t)
            k = sum(amps[i] * w[i] for i in range(6)) + (
                (tau - t).sin() if e == 0 else 0
            )
            v = expected[e] * k
            assert v > 0, (j, e, v)
            vals.append(v.lower())
        new = [w[0].lower(), w[1].lower(), min(vals)]
        mins = new if mins is None else [min(x, y) for x, y in zip(mins, new)]
    return list(map(str, mins))


def profile(p, cuts):
    Xs = [zero] + cuts + [pi]
    pt = (one / 2, p)
    data = []
    for i, (a, b) in enumerate(zip(Xs, Xs[1:])):
        r = int(i % 2 == 0)
        c = sub(pt, scale(unit(a), r))
        data.append((r, c))
        pt = add(c, scale(unit(b), r))

    def idx(t):
        z = [i for i, (a, b) in enumerate(zip(Xs, Xs[1:])) if a < t < b]
        assert len(z) == 1
        return z[0]

    def hvec(v):
        a = arb.atan2(v[1], v[0])
        assert a > 0 or a < 0
        sign = 1 if a > 0 else -1
        t = a if sign == 1 else a + pi
        r, c = data[idx(t)]
        return (r if sign == 1 else 1 - r) + dot(c, v)

    ev = sorted(
        [zero, L, *[x - (x / L).floor().unique_fmpz() * L for x in cuts]],
        key=lambda x: float(x.mid()),
    )
    rows = []
    assert all(a < b for a, b in zip(ev, ev[1:]))
    for i, (a, b) in enumerate(zip(ev, ev[1:])):
        c = -one / 2
        V = (zero, zero)
        for j in range(3):
            r, v = data[idx((a + b) / 2 + j * L)]
            s = (-1) ** j
            c += s * r
            V = add(V, scale(turn(v, -j * L), s))
        val = lambda y: c + dot(V, unit(y))
        der = lambda y: dot(V, perp(unit(y)))
        da, db = der(a), der(b)
        assert (da > 0 or da < 0) and (db > 0 or db < 0)
        if i == 0:
            assert da < 0 and val(a).contains(0)
        else:
            assert val(a) < 0
        if i == len(ev) - 2:
            assert db > 0 and val(b).contains(0)
        else:
            assert val(b) < 0
        peak = da > 0 and db < 0
        if peak and i != 1:
            assert c + magnitude(V) < 0
        if i == 1:
            assert peak
        rows.append(
            {
                "interior_maximum": peak,
                "left_F_prime": str(da),
                "right_F_prime": str(db),
            }
        )
    others = []
    sigma = rat(Q("1.5497")) * pi / 180
    for k in range(0, 360, 60):
        for sign in (1, -1):
            if (k, sign) in ((0, 1), (240, -1)):
                continue
            choice = 3 if SELECT[len(others)] == 0 else 5
            v = unit(pi / 6 + sigma + choice * L)
            u = turn((v[0], sign * v[1]), -k * pi / 180)
            gap = hvec(u) - one / 2
            assert gap > 0
            others.append(gap)
    return rows, str(min(x.lower() for x in others))


bands = [(Q(".4145"), Q(".4157"), 1), (Q(".442"), Q(".4442"), -1)]


def tangent_root(A, kind):
    Y = G.qa(A)
    brackets = {
        "peak": [
            ("0.41567220950", "0.41567220956"),
            ("0.44202986184", "0.44202986189"),
        ],
        "fixed": [
            ("0.41559931513", "0.41559931520"),
            ("0.44208109388", "0.44208109395"),
        ],
    }[kind]
    records = []
    for (l, r, s), (ll, rr) in zip(bands, brackets):
        threshold = [Q(ll), Q(rr)]
        assert l < threshold[0] < threshold[1] < r
        ends = []
        for w in threshold:
            q, _, u, _ = jets(G, rat(w))
            ends.append(dot(sub(Y, q), u))
        assert s * ends[0] < 0 < s * ends[1], (ll, rr, ends)
        worst = None
        N = 512
        for i in range(N):
            lb = l + (r - l) * i / N
            rb = l + (r - l) * (i + 1) / N
            mid = (lb + rb) / 2
            radius = rat((rb - lb) / 2)
            x = rat(lb).union(rat(rb))
            q, v, u, du = jets(G, rat(mid))
            qb, vb, ub, dub = jets(G, x)
            # Differentiate once more for a narrow rigorous derivative using q' perpendicular to u.
            w = arb_series([x, 1, 0, 0], prec=4)
            a, b = G.ns["points"](w)[0][2]
            xy = G.qa((a, b))
            assert all(z.prec >= 4 for z in xy)
            d = [z.derivative() for z in xy]
            A2 = d[0] ** 2 + d[1] ** 2
            speed = A2.sqrt()
            unitseries = [d[1] / speed, -d[0] / speed]
            assert all(z.prec >= 3 for z in unitseries)
            udd = tuple(2 * z[2] for z in unitseries)
            second = norm_up(vb) * norm_up(dub) + norm_up(sub(Y, qb)) * norm_up(udd)
            deriv = (s * dot(sub(Y, q), du) - radius * second).lower()
            assert deriv > 0, (i, deriv)
            direction = G.qa(u)
            angle = arb.atan2(direction[1], direction[0])
            angle_error = norm_up(dub) * radius
            assert (angle - angle_error) > rat(Q("4.8")) * pi / 180 and (
                angle + angle_error
            ) < rat(Q("5.1")) * pi / 180
            worst = deriv if worst is None else min(worst, deriv)
        records.append(
            {
                "threshold_bracket": [str(ll), str(rr)],
                "sign_at_bracket_ends": list(map(str, ends)),
                "derivative_orientation": s,
                "oriented_derivative_lower": str(worst),
                "certified_cells": N,
            }
        )
    return records


out = {}
for kind, pc, xc in [
    ("fixed", "0.13913108097436", "626.59561749555"),
    ("peak", "0.13913110314514", "626.59435541715"),
]:
    p, X, eu, A, V, theta, cuts, contraction = critical(kind, rat(Q(pc)), rat(Q(xc)))
    wc = calibrate(theta, cuts)
    profiles, other = profile(p, cuts)
    rad = dot(V, V) - one / 4
    if kind == "peak":
        assert rad.contains(0)
    else:
        assert rad > 0
    out[kind] = {
        "p": str(p),
        "x2_in_bin_units": str(X),
        "initial_contraction": str(contraction),
        "theta_degrees": str(theta * 180 / pi),
        "event_kernel_and_multiplier_lower": wc,
        "paired_margin_roots": tangent_root(A, kind),
        "F_profile_cells": profiles,
        "other_ten_cap_violation_lower": other,
        "critical_F_peak": str(magnitude(V) - one / 2),
    }
res = {
    "variants": out,
    "python": platform.python_version(),
    "python_flint": flint.__version__,
    "precision_bits": ctx.prec,
    "certifier_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
    "upstream_geometry_hashes": PINS,
    "dual_analytic_sha256": sha256((src / "dual_analytic.py").read_bytes()).hexdigest(),
}
args.output.write_text(json.dumps(res, indent=2, sort_keys=True) + "\n")
print(
    [
        (
            k,
            v["theta_degrees"],
            v["paired_margin_roots"][0]["threshold_bracket"],
            v["paired_margin_roots"][1]["threshold_bracket"],
        )
        for k, v in out.items()
    ]
)
