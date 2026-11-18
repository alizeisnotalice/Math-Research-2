"""Supporting fan velocities in nine directions transverse to the five-body manifold."""

from argparse import ArgumentParser
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import platform

import flint
from flint import arb, arb_mat, ctx
from certify_five_body_contact_stationary import (
    BOX,
    add,
    model,
    det,
    dot,
    interval,
    left,
    normal,
    q,
    r,
    scale,
    sub,
)

EVENT_ORDER = ("A", "Aopp", "B", "Bopp", "D", "F")


def plus(a, b):
    return [x + y for x, y in zip(a, b)]


def minus(a, b):
    return [x - y for x, y in zip(a, b)]


def rigid_velocities(owner, point, t):
    v = [(0, 0) for _ in range(9)]
    if owner == 5:
        v[:3] = [left(sub(point, t[3:5])), (1, 0), (0, 1)]
    elif owner == 7:
        v[3:6] = [left(sub(point, t[6:8])), (1, 0), (0, 1)]
    elif owner == 9:
        v[6:9] = [left(sub(point, t[9:11])), (1, 0), (0, 1)]
    return v


def connector(i, P, j, Q, t):
    W = map(add, rigid_velocities(i, P, t), rigid_velocities(j, Q, t))
    return [det(v, sub(Q, P)) / 2 for v in W]


def face_orders(geometry, ns):
    p = geometry
    a, b, f, d = (ns[k] for k in ("A", "B", "F", "D"))
    M = {
        "A": [
            (0, scale(r, a)),
            (3, p[3][0]),
            (5, add(p[5][3], a)),
            (7, p[7][0]),
            (9, p[9][0]),
        ],
        "B": [(0, scale(r, b)), (3, p[3][1]), (5, p[5][1]), (7, p[7][1]), (9, p[9][2])],
        "F": [(0, scale(r, f)), (5, add(p[5][3], f)), (9, p[9][1])],
        "D": [(0, scale(r, d)), (3, p[3][0]), (7, add(p[7][2], d)), (9, p[9][8])],
    }
    for label, points in M.items():
        n = ns[label]
        tangent = left(n)
        points.sort(key=lambda row: float(dot(tangent, row[1]).mid()))
        assert all(
            dot(tangent, p2) > dot(tangent, p1)
            for (_, p1), (_, p2) in zip(points, points[1:])
        )
    for label in ("A", "B"):
        M[label + "opp"] = [(i, sub(v, ns[label])) for i, v in reversed(M[label])]
    return M


def gradients(a, b):
    pieces, t, p, _ = model(a, b)
    ns = {tag: normal(stop) for tag, stop, *_ in pieces if tag in EVENT_ORDER}
    orders = face_orders(p, ns)
    base = [arb(0)] * 9
    increments = {}
    start = arb(0)
    for (tag, stop, i, _, rad, center), (_, _, j, _, rad2, center2) in zip(
        pieces, pieces[1:] + pieces[:1]
    ):
        n0, n1 = normal(start), normal(stop)
        delta = sub(n1, n0)
        base = plus(base, [rad * det(v, delta) for v in rigid_velocities(i, center, t)])
        P = add(center, scale(rad, n1))
        Q = add(center2, scale(rad2, n1))
        g = connector(i, P, j, Q, t)
        base = plus(base, g)
        if tag in orders:
            seq = orders[tag]
            assert seq[0][0] == i and seq[-1][0] == j
            increments[tag] = []
            for select in product((0, 1), repeat=len(seq) - 2):
                items = [
                    seq[0],
                    *(item for yes, item in zip(select, seq[1:-1]) if yes),
                    seq[-1],
                ]
                inc = [arb(0)] * 9
                for old, new in zip(items, items[1:]):
                    inc = plus(inc, connector(*old, *new, t))
                increments[tag].append(minus(inc, g))
        start = stop
    return (
        base,
        [increments[k] for k in EVENT_ORDER],
        {k: [i for i, _ in orders[k]] for k in EVENT_ORDER},
    )


def locate():
    import numpy as np
    from scipy.optimize import linprog

    base, inc, order = gradients(*(interval(*row) for row in BOX))
    print("order", order)
    b = np.asarray(list(map(float, base)))
    v = [np.asarray([list(map(float, r)) for r in ev]) for ev in inc]
    choices = np.array(list(product(*(range(len(ev)) for ev in v))))
    gs = b + sum(ev[choices[:, k]] for k, ev in enumerate(v))
    rng = np.random.default_rng(20261118)
    best = None
    for trial in range(400):
        ix = rng.choice(len(gs), 768, replace=False)
        sub = gs[ix]
        A = np.r_[np.ones((1, len(sub))), sub.T]
        rhs = np.r_[1.0, np.zeros(9)]
        ans = linprog(
            rng.normal(size=len(sub)),
            A_eq=A,
            b_eq=rhs,
            bounds=(0, None),
            method="highs",
        )
        if not ans.success:
            continue
        ids = np.flatnonzero(ans.x > 1e-9)
        if len(ids) != 10:
            continue
        Bi = np.linalg.inv(A[:, ids])
        w = Bi[:, 0]
        margin = min(w / np.sum(abs(Bi[:, 1:]), axis=1))
        if best is None or margin > best[0]:
            best = (margin, choices[ix[ids]].tolist(), min(w))
    print("best", best)


SIMPLEX = (
    (4, 3, 2, 5, 1, 0),
    (7, 5, 1, 7, 0, 1),
    (3, 1, 5, 3, 2, 0),
    (7, 5, 2, 1, 3, 0),
    (6, 5, 2, 6, 1, 1),
    (5, 4, 5, 2, 1, 0),
    (2, 1, 1, 1, 1, 0),
    (0, 5, 6, 6, 0, 0),
    (2, 1, 0, 6, 0, 0),
    (3, 5, 4, 0, 0, 0),
)


def certify():
    nodes = tuple(interval(*row) for row in BOX)
    pieces, t, p, ro = model(*nodes)
    A, B = normal(nodes[0]), normal(nodes[1])
    assert det(A, B) > 0 and all(
        abs(dot(left(v), p[3][k])) > q(".1") for k, v in enumerate([A, B])
    )
    base, incs, orders = gradients(*nodes)
    fans = []
    for cols in SIMPLEX:
        g = list(base)
        for j, k in enumerate(cols):
            g = plus(g, incs[j][k])
        fans.append(g)
    M = arb_mat([[1] * 10, *[[g[k] for g in fans] for k in range(9)]])
    assert M.det() > 0 or M.det() < 0
    inv = M.inv()
    eps = q(".00027")
    for i in range(10):
        assert inv[i, 0] > q(".05") and inv[i, 0] > eps * sum(
            abs(inv[i, j]) for j in range(1, 10)
        )
    path = Path(__file__)
    print(
        json.dumps(
            dict(
                precision_bits=ctx.prec,
                python=platform.python_version(),
                flint=flint.__version__,
                source_sha256=sha256(path.read_bytes()).hexdigest(),
                geometry_sha256=sha256(
                    path.with_name(
                        "certify_five_body_contact_stationary.py"
                    ).read_bytes()
                ).hexdigest(),
                normal_box=BOX,
                sorted_orders=orders,
                available_fans=32768,
                selected_choices=SIMPLEX,
                simplex_determinant=str(M.det()),
                barycentric_coordinates=[str(inv[i, 0]) for i in range(10)],
                transverse_gradient_cube_radius="27/100000",
                scope="Strict full eleven-dimensional local minimum at the isolated contact point",
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--locate", action="store_true")
    locate() if parser.parse_args().locate else certify()
