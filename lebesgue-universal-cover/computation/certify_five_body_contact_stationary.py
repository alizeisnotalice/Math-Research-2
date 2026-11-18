"""Two-angle contact chart for the five-body opt11 placement."""

import json, platform
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
import flint
from flint import arb, arb_series, ctx
from certify_four_body_contact_stationary import (
    add,
    sub,
    dot,
    det,
    scale,
    left,
    normal,
    lin2,
    unit_root,
    other,
    regular,
    r,
    pi,
    q,
    interval,
    primitive_key,
    ratio_angle,
)

ctx.prec = 192
ctx.cap = 3


def pose(n, A, B, N, spec, sign):
    v, w, z = (regular(n, i) for i, h in spec)
    (i, ha), (j, hb), (k, hn) = spec
    off = lin2(A, B, ha, hb)
    tc = lin2(A, B, -dot(A, v), -dot(B, w))
    ts = lin2(A, B, -dot(A, left(v)), -dot(B, left(w)))
    coeff = (dot(N, add(z, tc)), dot(N, add(left(z), ts)))
    cs, sn = unit_root(coeff, hn - dot(N, off), sign)
    if isinstance(cs, arb):
        assert cs > 0
    trans = add(off, add(scale(cs, tc), scale(sn, ts)))
    corners = []
    for k in range(n):
        v = regular(n, k)
        corners.append(add(trans, add(scale(cs, v), scale(sn, left(v)))))
    return (sn / cs).atan(), trans, corners


def model(a, b):
    A, B = normal(a), normal(b)
    t3 = lin2(A, B, r - dot(A, regular(3, 0)), r - dot(B, regular(3, 1)))
    p3 = [add(t3, regular(3, i)) for i in range(3)]
    D = other(p3[0], r, A)
    f7, t7, p7 = pose(7, A, B, D, [(0, r), (1, r), (2, -r)], -1)
    f9, t9, p9 = pose(9, A, B, D, [(0, r), (2, r), (8, r)], 1)
    F = unit_root(p9[1], r, -1)
    f5, t5, p5 = pose(5, A, B, F, [(3, -r), (1, r), (3, -r)], 1)
    p = {3: p3, 5: p5, 7: p7, 9: p9}
    us = {
        "C": other(p7[2], -r, D),
        "D": D,
        "E": other(p7[0], r, A),
        "F": F,
        "X": unit_root(sub(p5[2], p3[0]), 1, -1),
        "Y": unit_root(p5[2], r, 1),
        "Z": other(p9[8], -r, scale(-1, D)),
        "U": unit_root(sub(p9[4], p7[3]), 0, 1),
        "V": unit_root(sub(p5[2], p7[2]), 0, 1),
    }
    t = (*t3, f5, *t5, f7, *t7, f9, *t9)
    turns = {"C": 2, "D": 2, "E": 0, "F": 0, "X": 1, "Y": 1, "Z": 1, "U": 1, "V": 2}
    ev = {k: ratio_angle(v) + turns[k] * pi for k, v in us.items()}
    ev.update(A=a, B=b, Aopp=pi + a, Bopp=pi + b, end=2 * pi)
    signs = {k: 1 if v % 2 == 0 else -1 for k, v in turns.items()}
    if isinstance(a, arb):
        assert all(signs[k] * us[k][0] > 0 for k in us)
    rot = {3: 0 * a, 5: f5, 7: f7, 9: f9}
    bounds = []
    for tag, n, k, rad in WORD:
        if isinstance(tag, tuple):
            u, j, s = tag
            end = rot[u] + 2 * pi * j / u + s * pi / (2 * u)
        else:
            end = ev[tag]
        v = p[n][k] if n else (0 * a, 0 * a)
        bounds.append((tag, end, n, k, rad, v))
    return bounds, t, p, rot


WORD = [
    ("A", 3, 0, 0),
    ("E", 7, 0, 0),
    ("F", 0, 0, r),
    ((5, 1, -1), 5, 3, 1),
    ("B", 5, 1, 0),
    ((3, 1, 1), 3, 1, 0),
    ("X", 3, 0, 1),
    ("Y", 5, 2, 0),
    ("Z", 0, 0, r),
    ((9, 4, -1), 9, 8, 1),
    ("U", 9, 4, 0),
    ((7, 3, 1), 7, 3, 0),
    ("Aopp", 7, 0, 1),
    ((3, 2, -1), 3, 0, 1),
    ((3, 2, 1), 3, 2, 0),
    ("Bopp", 3, 1, 1),
    ((5, 4, -1), 5, 1, 1),
    ((5, 4, 1), 5, 4, 0),
    ("V", 5, 2, 1),
    ("C", 7, 2, 1),
    ("D", 0, 0, r),
    ("end", 3, 0, 0),
]


def area(a, b):
    pieces, *_ = model(a, b)
    start = 0 * a
    value = 0 * a
    for row, nxt in zip(pieces, pieces[1:] + pieces[:1]):
        tag, end, n, k, rad, v = row
        rad2 = nxt[4]
        w = nxt[5]
        n0, n1 = normal(start), normal(end)
        value += rad * det(v, sub(n1, n0)) + rad**2 * (end - start)
        value += det(add(v, scale(rad, n1)), add(w, scale(rad2, n1)))
        start = end
    return value / 2


def jets(a, b):
    x = area(arb_series([a, 1]), arb_series([b]))
    y = area(arb_series([a]), arb_series([b, 1]))
    z = area(arb_series([a, 1]), arb_series([b, 1]))
    h = z[2] - x[2] - y[2]
    return (x[1], y[1]), ((2 * x[2], h), (h, 2 * y[2]))


TIES = {k: {0, 3, 5, 7, 9} for k in ("A", "B", "Aopp", "Bopp")}
TIES.update(
    {
        "D": {0, 3, 7, 9},
        "F": {0, 5, 9},
        "C": {0, 7},
        "E": {0, 7},
        "X": {3, 5},
        "Y": {0, 5},
        "Z": {0, 9},
        "U": {7, 9},
        "V": {5, 7},
    }
)
BOX = (
    (".51344083304023783470819109772128672", ".51344083304023783470819109772128674"),
    ("1.57829241385673393299924195305530265", "1.57829241385673393299924195305530268"),
)


def verify_topology(a, b):
    pieces, t, p, rotations = model(a, b)
    cuts = {"start": arb(0), "end": 2 * pi}
    for n, f in rotations.items():
        for i in range(n):
            for end in (-1, 1):
                angle = f + 2 * pi * i / n + end * pi / (2 * n)
                if angle < 0:
                    angle += 2 * pi
                if angle > 2 * pi:
                    angle -= 2 * pi
                assert 0 < angle < 2 * pi
                cuts[n, i, end] = angle
    for k, angle, *_ in pieces:
        if k not in cuts:
            cuts[k] = angle
    keys = sorted(cuts, key=lambda k: float(cuts[k].mid()))
    assert [*keys[:1], *keys[-1:]] == ["start", "end"]
    assert all(cuts[k] < cuts[j] for k, j in zip(keys, keys[1:]))
    active = 0
    endpoints = 0
    critical = 0
    for k, j in zip(keys, keys[1:]):
        l, h = cuts[k], cuts[j]
        assert h - l < pi
        data = pieces[active][2:5]
        midpoint = (l + h) / 2
        if data[0] in rotations:
            assert primitive_key(data[0], midpoint, rotations[data[0]]) == data
        owner, vertex, rad = data
        center = p[owner][vertex] if owner else (arb(0), arb(0))
        nl, nh = normal(l), normal(h)
        for n in (0, 3, 5, 7, 9):
            if owner == n:
                continue
            if n:
                _, j2, r2 = primitive_key(n, midpoint, rotations[n])
                c2 = p[n][j2]
            else:
                r2 = r
                c2 = (arb(0), arb(0))
            v = sub(center, c2)
            equality = frozenset((n, owner))
            for tag, ll in [(k, nl), (j, nh)]:
                if equality <= TIES.get(tag, set()):
                    endpoints += 1
                else:
                    gap = rad - r2 + dot(v, ll)
                    assert gap > 0, (tag, owner, n, gap)
            # An interior minimum can occur only at normal -v/|v|.
            nv = scale(-1, v)
            if not (det(nl, nv) < 0 or det(nv, nh) < 0):
                margin = rad - r2 - dot(v, v).sqrt()
                assert margin > 0, (k, j, owner, n, margin)
                critical += 1
        if j == pieces[active][0]:
            active += 1
    assert active == len(pieces)
    return {
        "atomic_normal_intervals": len(keys) - 1,
        "exact_contact_endpoint_comparisons": endpoints,
        "interior_minimum_comparisons": critical,
    }


OLD = [
    "-0.0021161707485670258",
    "-0.002171588117597346",
    "0.03368534128826296",
    "0.025084694928075076",
    "-0.0038913081911800527",
    "0.708816831902889",
    "0.0035625741962694224",
    "-0.012629515020220292",
    "0.3490932487151976",
    "-0.0006351709996987754",
    "-0.0006510120748306202",
]


def improve_previous(new_area):
    # Decimal entries of independent/opt11_point.json, treated as exact rationals.
    previous = list(map(q, OLD))
    A = unit_root(add(previous[:2], regular(3, 0)), r, 1)
    B = unit_root(add(previous[:2], regular(3, 1)), r, -1)
    assert A[0] > 0 and B[0] < 0
    a, b = ratio_angle(A), pi + ratio_angle(B)
    verify_topology(a, b)
    anchor = model(a, b)[1]
    ds = [abs(z - w) for z, w in zip(anchor, previous)]
    # Rotation changes a support value by at most tan(pi/(2n))*|theta|/2.
    delta = ds[0] + ds[1]
    for n, i in [(5, 2), (7, 5), (9, 8)]:
        d = ds[i] * (pi / (2 * n)).tan() / 2 + ds[i + 1] + ds[i + 2]
        delta = arb.max(delta, d)
    # Bound the old perimeter by the circumference of its containing disk.
    R = q(".62")
    for n, tx, ty in [
        (3, *previous[:2]),
        (5, *previous[3:5]),
        (7, *previous[6:8]),
        (9, *previous[9:11]),
    ]:
        assert abs(tx) + abs(ty) + 1 / (2 * (pi / (2 * n)).cos()) < R
    old_lower = area(a, b) - 2 * pi * R * delta - pi * delta**2
    gain = old_lower - new_area
    assert gain > q(".00000000012")
    return {
        "old_parameters": OLD,
        "old_source_sha256": "4a8138a5d6014894fb2cad73b40fe3430a5257054558e92a7646713d8761b697",
        "hausdorff_upper": str(delta),
        "old_area_lower": str(old_lower),
        "guaranteed_area_decrease": "0.00000000012",
    }


def certify():
    X = [interval(*row) for row in BOX]
    mid = [q((F(row[0]) + F(row[1])) / 2) for row in BOX]
    rad = [q((F(row[1]) - F(row[0])) / 2) for row in BOX]
    g, H = jets(*X)
    g0, H0 = jets(*mid)
    assert H[0][0] > 0 and det(*H) > 0
    d = det(*H0)
    inv = ((H0[1][1] / d, -H0[0][1] / d), (-H0[1][0] / d, H0[0][0] / d))
    C0 = tuple(
        tuple(arb(round(float(v.mid()) * 10**12)) / 10**12 for v in row) for row in inv
    )
    assert det(*C0) > 0 or det(*C0) < 0
    defect = [
        [int(i == j) - sum(C0[i][k] * H[k][j] for k in range(2)) for j in range(2)]
        for i in range(2)
    ]
    contraction = max(sum(abs(v).abs_upper() for v in row) for row in defect)
    assert contraction < q(".001")
    for i in range(2):
        deviation = abs(dot(C0[i], g0)) + sum(
            abs(v) * rr for v, rr in zip(defect[i], rad)
        )
        assert deviation < rad[i]
    topo = verify_topology(*X)
    value = area(*mid)
    allowance = sum(abs(dg) * rr for dg, rr in zip(g0, rad)) + sum(
        abs(H[i][j]) * rad[i] * rad[j] / 2 for i in range(2) for j in range(2)
    )
    value = arb(value, allowance.abs_upper())
    assert q(".83729500622102524") < value < q(".83729500622102525")
    pieces, t, _, _ = model(*X)
    print(
        json.dumps(
            {
                "python": platform.python_version(),
                "flint": flint.__version__,
                "precision_bits": ctx.prec,
                "script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
                "stationary_angle_box": BOX,
                "contraction_norm_upper": str(contraction),
                "Hessian": [list(map(str, row)) for row in H],
                "Hessian_determinant": str(det(*H)),
                "placement": list(map(str, t)),
                "support_contact_angles_rad": {
                    str(k): str(ang) for k, ang, *_ in pieces if k in TIES
                },
                "verified_topology": topo,
                "improvement_over_opt11": improve_previous(value),
                "area_ball": str(value),
                "strict_area_bracket": [".83729500622102524", ".83729500622102525"],
                "scope": "Stationary minimum on the verified two-angle contact manifold, full placement admissible.",
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    certify()
