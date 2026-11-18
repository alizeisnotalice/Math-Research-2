"""An all-angle translation-dual certificate for the fixed seven-body placement."""

from argparse import ArgumentParser
from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
from math import ceil
from pathlib import Path
import json
import numpy as np
from scipy.optimize import linprog
from scipy.spatial import ConvexHull
from flint import arb, ctx
from harmonic_inner_material import polygon, det, ar, S, DEN, PI_QUALITY

ctx.prec = PI_QUALITY
SOURCE = Path(__file__).with_name("056-free-seven.json")
TGRID = 10**9


@lru_cache(None)
def directions(m):
    out = []
    for j in range(m):
        th = 2 * j * arb.pi() / m
        s, c = th.sin_cos()
        s11, c11 = (11 * th).sin_cos()
        out.append((s, c, s11, c11))
    return out


def material(index, count, m):
    phase = (2 * index + 1) * arb.pi() / count
    sp, cp = phase.sin_cos()
    # Center-phase and true-phase supports differ by <= (1/240)*pi/count.
    diff = Q(22, 7 * 240 * count)
    shrink = S - ceil(S * (diff + Q(2, S)) / Q(49, 100))
    pts = []
    for s, c, s11, c11 in directions(m):
        h = ar(Q(1, 2)) + (c11 * cp + s11 * sp) / 240
        hp = -(s11 * cp - c11 * sp) * 11 / 240
        X, Y = (h * c - hp * s) * S, (h * s + hp * c) * S
        x, y = round(float(X.mid())), round(float(Y.mid()))
        assert x - 1 < X < x + 1 and y - 1 < Y < y + 1
        pts.append((x * shrink, y * shrink))
    return pts


def numeric_plane(P, K, site):
    v = site / TGRID
    pts = np.vstack((P, K + v))
    cv = ConvexHull(pts)
    ix = cv.vertices
    y = pts[ix]
    d = np.roll(y, -1, axis=0) - np.roll(y, 1, axis=0)
    d = d[ix >= len(P)].sum(axis=0) / 2
    g = np.array([d[1], -d[0]])
    return cv.volume - g @ v, g


def exact_plane(P, K, site):
    tx, ty = (int(k) * (DEN // TGRID) for k in site)
    pts = sorted([(x, y, 0) for x, y in P] + [(x + tx, y + ty, 1) for x, y in K])
    assert all(p[:2] != q[:2] for p, q in zip(pts, pts[1:]))

    def turn(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    def half(order):
        out = []
        for p in order:
            while len(out) > 1 and turn(out[-2], out[-1], p) < 0:
                out.pop()
            out.append(p)
        return out[:-1]

    q = half(pts) + half(reversed(pts))
    # Collinear boundary points of one source body can be dropped without
    # losing differentiability. A mixed collinearity would need separate handling.
    for i, b in enumerate(q):
        a, c = q[i - 1], q[(i + 1) % len(q)]
        assert turn(a, b, c) > 0 or a[2] == b[2] == c[2]
    twice = sum(det(q[i - 1], q[i]) for i in range(len(q)))
    gx = sum(q[(i + 1) % len(q)][1] - q[i - 1][1] for i in range(len(q)) if q[i][2])
    gy = sum(q[i - 1][0] - q[(i + 1) % len(q)][0] for i in range(len(q)) if q[i][2])
    return twice - gx * tx - gy * ty, (gx, gy)


def lower(planes):
    (v0, g0), (v1, g1), (v2, g2) = planes
    w = [det(g1, g2), det(g2, g0), det(g0, g1)]
    if max(w) <= 0:
        w = [-z for z in w]
    if min(w) < 0 or sum(w) == 0:
        return Q(0)
    assert sum(a * b[0] for a, b in zip(w, [g0, g1, g2])) == 0
    assert sum(a * b[1] for a, b in zip(w, [g0, g1, g2])) == 0
    return Q(sum(a * b for a, b in zip(w, [v0, v1, v2])), sum(w) * 2 * DEN**2)


def suggest(P, K, target):
    p = np.array(P, dtype=float) / DEN
    k = np.array(K, dtype=float) / DEN
    sites = [
        np.array(z, dtype=int)
        for z in [
            (0, 0),
            (-30000000, -30000000),
            (-30000000, 30000000),
            (30000000, -30000000),
            (30000000, 30000000),
        ]
    ]
    planes = []
    for step in range(60):
        for s in sites[len(planes) :]:
            planes.append(numeric_plane(p, k, s))
        c = np.array([v for v, g in planes])
        g = np.array([g for v, g in planes])
        lp = linprog(
            [0, 0, 1],
            A_ub=np.c_[g, -np.ones(len(g))],
            b_ub=-c,
            bounds=[(None, None)] * 3,
            method="highs",
            options={
                "dual_feasibility_tolerance": 1e-10,
                "primal_feasibility_tolerance": 1e-10,
            },
        )
        assert lp.success
        if lp.fun > float(target) + 1e-11:
            ix = np.where(lp.ineqlin.marginals < -1e-10)[0]
            if len(ix) == 3:
                winner = [sites[i].tolist() for i in ix]
                if lower([exact_plane(P, K, z) for z in winner]) > target:
                    return winner
        s = np.rint(lp.x[:2] * TGRID).astype(int)
        if any(np.array_equal(s, t) for t in sites):
            raise RuntimeError(("uncertified", float(target), lp.fun))
        sites.append(s)
    raise RuntimeError(("iterations", lp.fun))


def main():
    ap = ArgumentParser()
    ap.add_argument("output", type=Path)
    ap.add_argument("--generate", action="store_true")
    ap.add_argument("--target", default=".83826")
    ap.add_argument("--quant", type=int, default=512)
    ap.add_argument("--angles", type=int, default=512)
    ap.add_argument("--degree", type=int, default=256)
    ap.add_argument("--test", type=int)
    a = ap.parse_args()
    raw = SOURCE.read_bytes()
    if a.generate:
        data = dict(
            source_sha256=sha256(raw).hexdigest(),
            target=a.target,
            angles=a.angles,
            quant=a.quant,
            base_degree=a.degree,
            translation_denominator=TGRID,
            rows=[],
        )
    else:
        data = json.loads(a.output.read_bytes())
    assert (
        data["source_sha256"] == sha256(raw).hexdigest()
        and data["translation_denominator"] == TGRID
    )
    target = Q(data["target"])
    P = polygon(raw, data["base_degree"])
    indices = [a.test] if a.test is not None else range(data["angles"])
    if not a.generate:
        assert len(data["rows"]) == data["angles"]
    worst = None
    for i in indices:
        K = material(i, data["angles"], data["quant"])
        if a.generate:
            row = suggest(P, K, target)
            data["rows"].append(row)
        else:
            row = data["rows"][i]
        lb = lower([exact_plane(P, K, s) for s in row])
        assert lb > target
        worst = lb if worst is None else min(lb, worst)
        if a.test is not None or i % 128 == 127:
            print(i, float(worst), flush=True)
    if a.generate and a.test is None:
        a.output.write_text(json.dumps(data, separators=(",", ":")) + "\n")
    print("verified", data["target"], len(indices), str(worst), flush=True)


if __name__ == "__main__":
    main()
