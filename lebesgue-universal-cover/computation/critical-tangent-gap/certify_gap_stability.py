"""Bound the Hausdorff loss of each calibration and place critical-band roots in the R212 core."""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import runpy
import sys

from flint import arb, arb_mat, ctx

parser = argparse.ArgumentParser()
parser.add_argument("computation", type=Path)
parser.add_argument("threshold_replay", type=Path)
parser.add_argument("output", type=Path)
args = parser.parse_args()
comp = args.computation.resolve()
base = Path(__file__).with_name("certify_gap_calibration.py")
# Record every upstream source and data file read by either proof, including transitive helpers.
read_sources = set()
local_files = {
    base.resolve(),
    Path(__file__).resolve(),
    args.output.resolve(),
    args.threshold_replay.resolve(),
}


def audit_read(event, data):
    if event != "open" or not isinstance(data[0], str):
        return
    path = Path(data[0])
    if path.suffix in (".py", ".json"):
        absolute = path.resolve()
        if absolute.is_relative_to(comp) and absolute not in local_files:
            read_sources.add(absolute.relative_to(comp))


sys.addaudithook(audit_read)
original = sys.argv[:]
try:
    sys.argv = [str(base), str(comp), str(args.threshold_replay)]
    ns = runpy.run_path(str(base), run_name="__main__")
finally:
    sys.argv = original
assert __debug__
rat, unit, dot, perp, G = (ns[k] for k in ("rat", "unit", "dot", "perp", "G"))
add, sub, scale, jets = (ns[k] for k in ("add", "sub", "scale", "jets"))
pi = arb.pi()
L = pi / 3
half = arb(1) / 2
deg = pi / 180
theta0 = rat(Q("6.1")) * deg
zero = arb(0)


def verify_stability(theta, cuts):
    cols = ns["sources"](theta)

    def row(t):
        return ns["analytic_row"](cols, t)

    pc = [sum((c * a.sin() for a, c in ev), zero) for ev, _, _ in cols]
    inverse = arb_mat([pc] + [row(x) for x in cuts]).inv()
    # Half-circle support query at alpha. Each modified kernel vanishes at all five switches.
    # Subdivide global alpha; on a straddling switch both the inactive and active branches are enclosed.
    bds = [zero] * 3
    angles = [zero, theta, L, theta + L, cols[0][0][0][0], 2 * L, theta + 2 * L, pi]
    assert all(a < b for a, b in zip(angles, angles[1:]))
    forms = []
    for a, b in zip(angles, angles[1:]):
        m = (a + b) / 2
        forms.append(
            (
                a,
                b,
                [
                    (
                        d1 - sum((c * x.cos() for x, c in ev if x > m), zero),
                        d2 + sum((c * x.sin() for x, c in ev if x > m), zero),
                    )
                    for ev, d1, d2 in cols
                ],
            )
        )
    for i in range(2048):
        t = pi * i / 2048
        alpha = t.union(pi * (i + 1) / 2048)
        rhs = [-alpha.sin()]
        for x in cuts:
            if alpha < x:
                rhs.append(zero)
            elif alpha > x:
                rhs.append(-(alpha - x).sin())
            else:
                rhs.append(zero.union(-(alpha - x).sin()))
        weight = inverse * arb_mat([[x] for x in rhs])
        limit = zero
        for a, b, form in forms:
            A = sum((c * weight[j, 0] for j, (c, d) in enumerate(form)), zero)
            B = sum((d * weight[j, 0] for j, (c, d) in enumerate(form)), zero)
            if not alpha > b:
                limit = max(limit, (A**2 + B**2).sqrt().upper())
            if not alpha < a:
                limit = max(
                    limit,
                    ((A - alpha.cos()) ** 2 + (B + alpha.sin()) ** 2).sqrt().upper(),
                )
        for j, value in enumerate((limit, abs(weight[0, 0]), abs(weight[1, 0]))):
            bds[j] = max(bds[j], value.upper())
    assert bds[0] < 6
    # On intervals where the master kernel has one sign, its absolute value is concave.
    # Voronoi midpoints make distance to the switches affine on each such interval.
    cosine = inverse * arb_mat([[0], [cuts[0].sin()], [0], [0], [0], [0]])
    sine = inverse * arb_mat([[-1], [-cuts[0].cos()], [0], [0], [0], [0]])
    low, high = rat(Q("4.8")) * deg, rat(Q("5.1")) * deg
    worst = None
    ah = cols[0][0][0][0]
    halfcuts = [(a + b) / 2 for a, b in zip(cuts, cuts[1:])]
    for i in range(256):
        a = low + (high - low) * i / 256
        tau = a.union(low + (high - low) * (i + 1) / 256)
        weights = [cosine[j, 0] * tau.cos() + sine[j, 0] * tau.sin() for j in range(6)]
        assert bds[1] < 1536 * weights[0] and bds[2] < 1536 * weights[1]
        for name, t in [
            ("0", zero),
            ("tau", tau),
            ("theta", theta),
            ("L", L),
            ("theta+L", theta + L),
            ("ah", ah),
            ("2L", 2 * L),
            ("theta+2L", theta + 2 * L),
            ("pi", pi),
        ] + [("bisector", y) for y in halfcuts]:
            if name == "tau":
                master = zero
            elif tau > t:
                master = (tau - t).sin()
            else:
                assert tau < t
                master = zero
            value = master + sum(k * v for k, v in zip(row(t), weights))
            assert all(x < t or x > t for x in cuts)
            sign = (-1) ** sum(x < t for x in cuts)
            dist = min(abs(t - x).upper() for x in cuts)
            quotient = (sign * value / dist).lower()
            assert quotient > arb(1) / 256, (name, quotient)
            worst = quotient if worst is None else min(worst, quotient)
    return dict(
        test_kernel_Lipschitz_upper=str(bds[0]),
        test_inequality_weights_upper=list(map(str, bds[1:])),
        master_kernel_distance_quotient_lower=str(worst),
        hausdorff_factor=1536,
    )


records = {}
for kind, pc, xc in [
    ("fixed", "0.13913108097436", "626.59561749555"),
    ("peak", "0.13913110314514", "626.59435541715"),
]:
    p, X, eu, A, V, theta, cuts, _ = ns["critical"](kind, rat(Q(pc)), rat(Q(xc)))
    records[kind] = verify_stability(theta, cuts)
    if kind == "fixed":
        p0, xs, A0 = p, [zero, *cuts, pi], A
    else:
        peak = theta


# Monotonicity from the calibration replay bounds the box under consideration by these outer endpoints.
def tangent_excess(w):
    q, _, u, _ = jets(G, rat(Q(w)))
    return dot(sub(G.qa(A0), q), u)


upper = max(
    tangent_excess("0.41567220956").upper(), tangent_excess("0.44202986184").upper()
)
eta = rat(Q("0.000009"))
assert 1536 * upper < eta
pieces = []
point = (half, p0)
for i, (a, b) in enumerate(zip(xs, xs[1:])):
    r = int(i % 2 == 0)
    c = sub(point, scale(unit(a), r))
    pieces.append((r, c))
    point = add(c, scale(unit(b), r))
core_inputs = {
    "audit_eh_continuous_placement.py": "25675cc2d431b80db5de491a9e649bb7a38941957c6b7214eb0179a7518b9dea",
    "audit_optimized_eh_chord.py": "dd371b0545a1e1ceef3edefba692623b592a904d4c16d372b715f9ffa6ca9c64",
    "eh-upper-chain-parameters.json": "053bdfc3781ec963155ad90807232b8f9d56bb3dd75cbe2f6b683ee41135a6d7",
    "certify_eh_upper_chain.py": "2748412ad86e0a421a515be67ae8dcf21dc89eafb7f51911b8a6525cda1773e3",
    "eh-two-line-cover-certificate.json": "10f2f3aff4eb8c60e869053bacf8f569d2db55dec305fc7c0b405cdecf3bddf2",
}
assert all(
    sha256((comp / name).read_bytes()).hexdigest() == digest
    for name, digest in core_inputs.items()
)
import audit_optimized_eh_chord as eh_geom
from audit_eh_continuous_placement import polygon
from certify_eh_upper_chain import front

ctx.prec = 224
eh_geom.SIGMA = rat(Q("1.5497")) * eh_geom.PI / 180
params = json.loads((comp / "eh-upper-chain-parameters.json").read_text())
old = json.loads((comp / "eh-two-line-cover-certificate.json").read_text())
eh = front(old, params)
assert len(eh["active"]) == 11
inner, _ = polygon(params[-1][0], 128)
hex_rows = {"B1-Q": 1, "P-F1": 2, "F1-H": 3, "E3-D1": 4, "D1-G": 5, "C3-B1": 0}
assert len(inner) == 906 and set(hex_rows) == {
    label for label, _, _ in inner if label in hex_rows
}
for j, (_, u, _) in enumerate(inner):
    _, v, _ = inner[(j + 1) % len(inner)]
    assert u[0] * v[1] - u[1] * v[0] > 0 and dot(u, v) > 0
for label, n, bound in inner:
    if label in hex_rows:
        v = unit(hex_rows[label] * L)
        assert bound.contains(arb(1) / 2) and all(x.overlaps(y) for x, y in zip(v, n))
ta = add(G.ns["z"], G.ns["D"])
tp = G.qa(G.ns["point_jets"](G.A(Q(1)))[0][2])
chord = sub(tp, ta)
core_normal = scale(perp(chord), -1 / dot(chord, chord).sqrt())
core_height = dot(core_normal, ta)
assert core_height > 0 and dot(core_normal, G.ns["z"]) > core_height

root3 = arb(3).sqrt()
R = 1 / root3
caps = [unit(pi / 6 + rat(Q("1.5497")) * deg + k * L) for k in (3, 5)]
rows = [("inner " + label, v, b) for label, v, b in inner if label not in hex_rows]
rows.extend(
    [
        ("EH " + str(i), v, b)
        for i, (v, b) in enumerate(zip(eh["normals"], eh["offsets"], strict=True))
    ]
)
rows.extend(
    [
        ("core chord", core_normal, core_height),
        ("AH", ns["G"].normal, ns["G"].height),
    ]
)


def coeffs(left, right, normal, phi=None):
    if phi is None:
        phi = arb.atan2(normal[1], normal[0])
    q = left.union(right) + phi
    nx, ny = normal
    forms = []
    for k in range(-2, 3):
        reduced = q - k * pi
        for i, (a, b) in enumerate(zip(xs, xs[1:])):
            if reduced < a or reduced > b:
                continue
            r, c = pieces[i]
            forms.append(
                (1 - r if k % 2 else r, c[0] * nx + c[1] * ny, c[1] * nx - c[0] * ny)
            )
    assert forms
    return forms


def block(left, right, normal):
    phi = arb.atan2(normal[1], normal[0])
    j = round(float((phi / L).mid()))
    gap = phi - j * L
    # Nearest face and the adjacent face on the side towards the row.
    s = 1 if float(gap.mid()) >= 0 else -1
    beta = gap.sin() / (s * L.sin())
    alpha = gap.cos() - beta * L.cos()
    v0 = unit(j * L)
    v1 = unit((j + s) * L)
    p0 = coeffs(left, right, v0, j * L)
    p1 = coeffs(left, right, v1, (j + s) * L)
    assert len(p0) == len(p1) == 1
    r0, a0, b0 = p0[0]
    r1, a1, b1 = p1[0]
    mid = (left + right) / 2
    rad = (right - left) / 2
    u = unit(mid)
    v = (-u[1], u[0])
    err = min(
        ((1 + abs(alpha) + abs(beta)) * eta).upper(),
        (
            ((2 * eta * (1 + eta)).sqrt()) * abs(gap)
            + (abs(1 - alpha) + abs(beta)) * eta
        ).upper(),
    )
    upper = []
    for r, a, b in coeffs(left, right, normal):
        D = r + alpha * (arb(1) / 2 - r0) + beta * (arb(1) / 2 - r1)
        A = a - alpha * a0 - beta * a1
        B = b - alpha * b0 - beta * b1
        bound = (
            D
            + A * u[0]
            + B * u[1]
            + abs(A * v[0] + B * v[1]) * rad
            + (abs(A) + abs(B)) * rad**2 / 2
            + err
        )
        upper.append(bound.upper())
    return max(upper)


rows.extend([("m3", caps[0], arb(1) / 2), ("m5", caps[1], arb(1) / 2)])
wlo, whi = Q(".405"), Q(".445")
N = 256
for k in range(N):
    a = wlo + (whi - wlo) * k / N
    b = wlo + (whi - wlo) * (k + 1) / N
    pos, _, u, _ = jets(G, rat((a + b) / 2))
    box, _, _, du = jets(G, rat(a).union(rat(b)))
    speed = sum(abs(x).upper() ** 2 for x in du).sqrt()
    normq = sum(abs(x).upper() ** 2 for x in box).sqrt()
    loss = rat((b - a) / 2) * (R + normq) * speed
    rows.append(("tangent " + str(k), u, dot(pos, u) - loss))

assert len(rows) == 1171
assert all(dot(n, n).contains(1) for _, n, _ in rows)


def support(t):
    inds = [i for i, (a, b) in enumerate(zip(xs, xs[1:])) if a < t < b]
    assert len(inds) == 1
    r, c = pieces[inds[0]]
    return r + dot(c, unit(t))


right = rat(Q("6.8")) * deg
Fright = support(right) - support(right + L) + support(right + 2 * L) - half
assert Fright < -3 * eta
fits = []
for left, end, Nangle in ((theta0, peak, 1), (peak, right, 12)):
    each = []
    for label, n, height in rows:
        ub = max(
            block(
                left + (end - left) * k / Nangle,
                left + (end - left) * (k + 1) / Nangle,
                n,
            )
            for k in range(Nangle)
        )
        each.append(((height - ub).lower(), label))
    sl, label = min(each)
    print(left / deg, end / deg, sl, label)
    assert sl > 0, (sl, label, sorted(each)[:4])
    fits.append(
        dict(
            degrees=list(map(str, (left / deg, end / deg))),
            angle_intervals=Nangle,
            minimum_label=label,
            minimum_slack_by_group={
                g: str(min(z for z, lab in each if lab.startswith(g)))
                for g in ("inner", "EH", "core chord", "AH", "m3", "m5", "tangent")
            },
        )
    )
out = dict(
    stability=records,
    fixed_tangent_budget_upper=str(upper),
    hausdorff_radius="0.000009",
    F_at_right_upper=str(Fright.upper()),
    pose_certificates=fits,
    core_dependency_hashes=core_inputs,
    threshold_replay_sha256=sha256(args.threshold_replay.read_bytes()).hexdigest(),
    script_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
    arb_precision_bits=ctx.prec,
)
for module in tuple(sys.modules.values()):
    origin = getattr(module, "__file__", None)
    if isinstance(origin, str):
        path = Path(origin).resolve()
        if (
            path.is_relative_to(comp)
            and path.suffix == ".py"
            and path not in local_files
        ):
            read_sources.add(path.relative_to(comp))
out["upstream_files_sha256"] = {
    path.as_posix(): sha256((comp / path).read_bytes()).hexdigest()
    for path in sorted(read_sources)
}
args.output.write_text(json.dumps(out, sort_keys=True, indent=2) + "\n")
