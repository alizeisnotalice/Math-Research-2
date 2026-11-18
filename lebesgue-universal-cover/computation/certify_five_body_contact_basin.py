"""Explicit curved and rectangular neighborhoods of the five-body strict minimum."""

from argparse import ArgumentParser
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import json
from flint import arb, arb_mat
from four_contact_taylor import Model
from certify_five_body_contact_stationary import (
    BOX,
    area,
    model,
    q,
    interval,
    normal,
    primitive_key,
    det,
    dot,
    left,
    add,
    sub,
    scale,
    pi,
    r,
    regular,
    unit_root,
)
from certify_five_body_contact_transverse import (
    EVENT_ORDER,
    SIMPLEX,
    face_orders,
    rigid_velocities,
    connector,
    plus,
)

CENTERS = (".513440833040237835", "1.578292413856733933")
RAW = (
    "-.002115486753128799",
    "-.002165784484548596",
    ".033675690138239897",
    ".025088017733157594",
    "-.003893364590431682",
    ".708851510983043453",
    ".003564914000574384",
    "-.012628775899428358",
    ".349093257175728252",
    "-.000634982923566646",
    "-.000649291654038299",
)


def evaluate(v, x, y):
    if isinstance(v, Model):
        return v.evaluate(x, y)
    if isinstance(v, (list, tuple)):
        return tuple(evaluate(k, x, y) for k in v)
    return v


def geometry_box(data, x, y):
    pcs, t, p, rots = data
    pcs = [
        (tag, evaluate(end, x, y), n, k, rho, evaluate(v, x, y))
        for tag, end, n, k, rho, v in pcs
    ]
    return (
        pcs,
        evaluate(t, x, y),
        {n: evaluate(v, x, y) for n, v in p.items()},
        {n: evaluate(v, x, y) for n, v in rots.items()},
    )


def atomic_order(pcs, rots):
    cuts = {"start": arb(0), "end": 2 * pi}
    for n, ang in rots.items():
        for k in range(n):
            for sg in (-1, 1):
                a = ang + 2 * pi * k / n + sg * pi / (2 * n)
                if a < 0:
                    a += 2 * pi
                if a > 2 * pi:
                    a -= 2 * pi
                assert 0 < a < 2 * pi
                cuts[n, k, sg] = a
    for key, a, *_ in pcs:
        cuts[key] = a
    order = sorted(cuts, key=lambda k: float(cuts[k].mid()))
    assert order[0] == "start" and order[-1] == "end"
    for k, j in zip(order, order[1:]):
        assert cuts[k] < cuts[j]
    return len(order) - 1


def links_and_arcs(data):
    pcs, t, p, rots = data
    ns = {tag: normal(a) for tag, a, *_ in pcs if tag in EVENT_ORDER}
    orders = face_orders(p, ns)
    arcs = []
    links = []
    prev = q(0)
    arc_grad = [arb(0)] * 9
    for row, nxt in zip(pcs, pcs[1:] + pcs[:1]):
        tag, end, i, k, rho, v = row
        _, _, j, _, rho2, w = nxt
        n0 = normal(prev)
        n1 = normal(end)
        diff = sub(n1, n0)
        if i in (5, 7, 9) and rho == 1:
            arcs.append((i, diff))
            arc_grad = plus(
                arc_grad, [det(vv, diff) for vv in rigid_velocities(i, v, t)]
            )
        P = add(v, scale(rho, n1))
        Q = add(w, scale(rho2, n1))
        if i != j:
            seq = orders.get(tag, [(i, P), (j, Q)])
            assert seq[0][0] == i and seq[-1][0] == j
            links.append((tag, seq))
        prev = end
    return t, links, arcs, arc_grad


def rotated(p, cosine, sine):
    return add(scale(cosine, p), scale(sine, left(p)))


def fan_bounds(t, links, arcs, d):
    delta = arb(0, d)
    co, si = delta.cos(), delta.sin()
    h = [[arb(0)] * 9 for _ in range(9)]
    e0, e1 = (1, 0), (0, 1)

    def moved(i, P):
        if i in (0, 3):
            return P, {}, {}
        k = 3 * ((i - 5) // 2)
        tr = t[k + 3 : k + 5]
        v = rotated(sub(P, tr), co, si)
        return (
            add(add(tr, (delta, delta)), v),
            {k: left(v), k + 1: e0, k + 2: e1},
            {k: scale(-1, v)},
        )

    for i, v in arcs:
        k = 3 * ((i - 5) // 2)
        trans = add(t[k + 3 : k + 5], (delta, delta))
        dv = rotated(v, co, si)
        h[k][k] -= det(trans, dv) / 2
        for j, e in enumerate((e0, e1), k + 1):
            h[k][j] += det(e, left(dv)) / 2
            h[j][k] = h[k][j]
    orientations = []
    for (i, P), (j, Q) in links:
        V, dv, hv = moved(i, P)
        W, dw, hw = moved(j, Q)
        orientations.append(det(V, W).lower())
        for k, dp in hv.items():
            h[k][k] += det(dp, W) / 2
        for k, dq in hw.items():
            h[k][k] += det(V, dq) / 2
        for k, dp in dv.items():
            for l, dq in dw.items():
                v = det(dp, dq) / 2
                h[k][l] += v
                h[l][k] += v
    return min(orientations), max(sum(abs(z).abs_upper() for z in row) for row in h)


def test(data, volume, R, x, y, d, rho):
    raw = geometry_box(data, x, y)
    pcs, t, p, rots = raw
    atomic_order(pcs, rots)
    H = [
        [volume.second(i, j, x, y, R) - q(".1") * (i == j) for j in range(2)]
        for i in range(2)
    ]
    assert H[0][0] > 0 and det(*H) > 0, ("hessian", det(*H))
    t, faces, arcs, arc_grad = links_and_arcs(raw)
    G = []
    orient = []
    tails = []
    for ch in SIMPLEX:
        sel = dict(zip(EVENT_ORDER, ch))
        chain = []
        g = list(arc_grad)
        for tag, seq in faces:
            if tag in sel:
                mid = len(seq) - 2
                seq = [
                    seq[0],
                    *(
                        u
                        for k, u in enumerate(seq[1:-1])
                        if sel[tag] & (1 << (mid - 1 - k))
                    ),
                    seq[-1],
                ]
            for u, v in zip(seq, seq[1:]):
                g = plus(g, connector(*u, *v, t))
                chain.append((u, v))
        orient1, sup = fan_bounds(t, chain, arcs, d)
        orient.append(orient1)
        tails.append(sup * d / 2)
        G.append(g)
    M = arb_mat([[1] * 10, *[[row[i] for row in G] for i in range(9)]])
    # A midpoint preconditioner prevents enclosure blowup from interval elimination.
    B = arb_mat([[M[i, j].mid() for j in range(10)] for i in range(10)]).inv()
    pre = arb_mat([[v.mid() for v in row] for row in B.tolist()])
    assert abs(pre.det()) > 0
    inv = (pre * M).inv() * pre
    slack = min(
        (inv[i, 0] - rho * sum(abs(inv[i, j]) for j in range(1, 10))).lower()
        for i in range(10)
    )
    linear = (rho - max(tails)).lower()
    min_o = min(orient)
    assert slack > 0 and linear > 0 and min_o > 0, ("transverse", slack, linear, min_o)
    # Every moved body contains the origin by its uniform inball and the small
    # translated-center norm. Thus the arcs are positively polar oriented.
    for i, n in enumerate((3, 5, 7, 9)):
        txy = t[3 * i : 3 * i + 2]
        assert abs(txy[0]) + abs(txy[1]) + 2 * d < 1 - 1 / (2 * (pi / (2 * n)).cos())
    return [det(*H).lower(), slack, min_o, linear]


def main():
    ap = ArgumentParser()
    ap.add_argument("--scale", default="0.000006")
    ap.add_argument("--transverse", default="0.00008")
    ap.add_argument("--cube", default="0.00012")
    ap.add_argument("--degree", type=int, default=9)
    ap.add_argument("--cells", type=int, default=20)
    ap.add_argument("--raw-triangle", default="0.00000025")
    ap.add_argument("--raw-other", default="0.000055")
    args = ap.parse_args()
    cx, cy = map(q, CENTERS)
    R = q(args.scale)
    d = q(args.transverse)
    Model.degree = args.degree
    a = Model({(0, 0): cx, (1, 0): R})
    b = Model({(0, 0): cy, (0, 1): R})
    data = model(a, b)
    volume = area(a, b)
    cells = args.cells
    assert (
        cells > 0
        and R > 0
        and d > 0
        and all(abs(interval(*row) - c) < R / 2 for row, c in zip(BOX, (cx, cy)))
    )
    margins = []
    for i in range(cells):
        for j in range(cells):
            unit = lambda k: arb(
                q(Fraction(2 * k + 1 - cells, 2 * cells)), q(Fraction(1, 2 * cells))
            )
            margins.append(test(data, volume, R, unit(i), unit(j), d, q(args.cube)))
    widths = [q(args.raw_triangle)] * 2 + [q(args.raw_other)] * 9
    raw = [arb(q(v), w) for v, w in zip(RAW, widths)]
    true_pose = model(*(interval(*row) for row in BOX))[1]
    assert all(
        abs(p - q(c)) < w for p, c, w in zip(true_pose, RAW, widths, strict=True)
    )
    z = raw[:2]
    A = unit_root(add(z, regular(3, 0)), r, 1)
    B = unit_root(add(z, regular(3, 1)), r, -1)
    assert A[0] > 0 and B[1] > 0
    x = ((A[1] / A[0]).atan() - cx) / R
    y = (pi / 2 - (B[0] / B[1]).atan() - cy) / R
    assert abs(x) < q(".5") and abs(y) < q(".5"), ("angles", x, y)
    anchor = evaluate(data[1], x, y)
    offsets = [abs(raw[i] - anchor[i]).abs_upper() for i in range(2, 11)]
    assert max(offsets) < d, ("deviations", list(map(str, offsets)))
    print(
        json.dumps(
            dict(
                dependencies={
                    k: sha256(Path(__file__).with_name(k).read_bytes()).hexdigest()
                    for k in (
                        "four_contact_taylor.py",
                        "certify_five_body_contact_stationary.py",
                        "certify_five_body_contact_transverse.py",
                    )
                },
                source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
                center=CENTERS,
                normal_halfwidth=str(Fraction(args.scale) / 2),
                transverse=args.transverse,
                gradient_cube=args.cube,
                raw_center=RAW,
                raw_halfwidth=[args.raw_triangle] * 2 + [args.raw_other] * 9,
                raw_normalized_angles=[str(x), str(y)],
                maximum_raw_deviation=str(max(offsets)),
                degree=args.degree,
                cells=cells**2,
                vol_error=str(volume.error),
                curvature_floor="0.1",
                interval_minima=[str(min(z[k] for z in margins)) for k in range(4)],
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
