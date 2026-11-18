"""Numerical proposal for a common placement of the seven bodies and an 11th harmonic."""

import json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.spatial import ConvexHull
from optimize_nine import shape_boundaries, unit


class Family:
    def __init__(self, k):
        bases = shape_boundaries(k, 8 * k, 16 * k)
        theta = np.arange(16 * k) * 2 * np.pi / (16 * k)
        u = unit(theta)
        h = 0.5 + np.cos(11 * theta) / 240
        hp = -11 * np.sin(11 * theta) / 240
        body = h[:, None] * u + hp[:, None] * np.c_[-u[:, 1], u[:, 0]]
        self.boundaries = bases[:5] + bases[7:] + [body]
        self.owners = np.concatenate(
            [np.full(len(v), i) for i, v in enumerate(self.boundaries)]
        )

    def value(self, x, derivative=True):
        points = [self.boundaries[0], self.boundaries[1] + x[:2]]
        for j, shape in enumerate(self.boundaries[2:]):
            r, tx, ty = x[2 + 3 * j : 5 + 3 * j]
            rot = np.array([[np.cos(r), np.sin(r)], [-np.sin(r), np.cos(r)]])
            points.append(shape @ rot + [tx, ty])
        p = np.vstack(points)
        H = ConvexHull(p)
        if not derivative:
            return H.volume
        q = p[H.vertices]
        owners = self.owners[H.vertices]
        dq = np.roll(q, -1, axis=0) - np.roll(q, 1, axis=0)
        g = np.c_[dq[:, 1], -dq[:, 0]] / 2
        xx, yy = [np.bincount(owners, weights=v, minlength=8) for v in g.T]
        angle = []
        for j in range(2, 8):
            idx = owners == j
            v = q[idx] - x[3 * (j - 2) + 3 : 3 * (j - 2) + 5]
            angle.append(sum(g[idx, 1] * v[:, 0] - g[idx, 0] * v[:, 1]))
        return H.volume, np.r_[
            xx[1], yy[1], np.stack([angle, xx[2:], yy[2:]], 1).ravel()
        ]


def main():
    from argparse import ArgumentParser

    ap = ArgumentParser()
    ap.add_argument("output", type=Path)
    args = ap.parse_args()
    x = np.asarray(
        json.loads(Path(__file__).with_name("056-free-seven.json").read_bytes())[
            "active_parameters"
        ],
        dtype=float,
    )
    f = Family(22)
    g = Family(100)
    hi = Family(200)
    sc = np.array([0.015, 0.015] + 6 * [0.09, 0.015, 0.015])
    seeds = []
    for theta in np.arange(10) * 2 * np.pi / 110:
        st = np.array([theta, 0.0, 0.0])
        a = minimize(
            lambda z: f.value(np.r_[x, z], False),
            st,
            method="Nelder-Mead",
            options=dict(
                maxiter=250,
                initial_simplex=np.vstack([st, st + np.diag([0.012, 0.003, 0.003])]),
            ),
        )
        seeds.append((g.value(np.r_[x, a.x], False), a.x.tolist()))
    allrows = []
    period = 2 * np.pi / 11
    distinct = []
    for _, st in sorted(seeds):
        if all(
            abs((st[0] - u[0] + period / 2) % period - period / 2) > 0.015
            for u in distinct
        ):
            distinct.append(st)
    for st in [u for _, u in sorted(seeds)[:4]] + distinct[1:4]:
        y = np.r_[x, st]
        for f0 in (f, g):

            def ev(z):
                area, grad = f0.value(y + sc * z)
                return area, grad * sc

            r = minimize(
                ev,
                sc * 0,
                method="L-BFGS-B",
                jac=True,
                bounds=[(-2, 2)] * len(sc),
                options=dict(maxiter=330, ftol=1e-15, gtol=1e-10, maxls=36),
            )
            y = y + sc * r.x
        a = hi.value(y, False)
        print("joint", a, flush=True)
        allrows.append(dict(sample_area=a, params=[format(k, ".17g") for k in y]))
        args.output.write_text(
            json.dumps(
                dict(
                    best=min(allrows, key=lambda r: r["sample_area"]), searched=allrows
                ),
                indent=2,
            )
            + "\n"
        )


if __name__ == "__main__":
    main()
