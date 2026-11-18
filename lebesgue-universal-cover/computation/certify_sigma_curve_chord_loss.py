"""Certified curvature separation and the inscribed polygon area loss on the compulsory AH paired arc."""

import argparse
import json
from pathlib import Path
from fractions import Fraction as F
from math import comb, factorial
from flint import arb, arb_series, ctx
from certify_sigma_mandatory_arcs import g, degree

ctx.cap = 10
n = [g.unit(i * g.PI / 3) for i in range(6)]
h = g.HALF
rt = arb(3).sqrt()
mid = (-h / 2, -rt / 4)
z = g.circle_pair((rt / 2 - h, (rt / 2 - h - 1) / rt), mid, -1)
el = (rt / 2, h)
en = (-rt / 2, h)


class SigmaJet:
    """Ordinary Taylor coefficients in sigma, each holding a w-series. Degree is at most three."""

    def __init__(self, v):
        self.v = v if isinstance(v, list) else [v]

    def __getitem__(self, i):
        return self.v[i] if i < len(self.v) else 0

    def __neg__(self):
        return SigmaJet([-x for x in self.v])

    def __add__(self, b):
        if not isinstance(b, SigmaJet):
            b = SigmaJet(b)
        return SigmaJet([self[i] + b[i] for i in range(max(len(self.v), len(b.v)))])

    __radd__ = __add__

    def __sub__(self, b):
        return self + (-b)

    def __rsub__(self, b):
        return b + (-self)

    def __mul__(self, b):
        if not isinstance(b, SigmaJet):
            return SigmaJet([x * b for x in self.v])
        ln = min(4, len(self.v) + len(b.v) - 1)
        return SigmaJet(
            [
                sum(
                    self[i] * b[k - i]
                    for i in range(max(0, k + 1 - len(b.v)), min(k + 1, len(self.v)))
                )
                for k in range(ln)
            ]
        )

    __rmul__ = __mul__

    def inverse(self):
        v = [1 / self[0]]
        for i in range(1, 4):
            v.append(-v[0] * sum(self[k] * v[i - k] for k in range(1, i + 1)))
        return SigmaJet(v)

    def __truediv__(self, b):
        if isinstance(b, SigmaJet):
            return self * b.inverse()
        return SigmaJet([x / b for x in self.v])

    def __rtruediv__(self, b):
        return self.inverse() * b

    def sqrt(self):
        v = [self[0].sqrt()]
        for i in range(1, 4):
            v.append((self[i] - sum(v[k] * v[i - k] for k in range(1, i))) / (2 * v[0]))
        return SigmaJet(v)

    def __pow__(self, k):
        assert k == 2
        return self * self


def slope(sig):
    """The direction T_a-Z as an interval Taylor series in sigma (radians)."""
    c = g.line_pair(n[4], g.unit(7 * g.PI / 6 + sig))
    c = (-c[0], c[1])
    v = g.minus(z, c)
    vv = g.dot(v, v)
    return g.minus(
        g.plus(c, g.scale(1 / vv, g.minus(v, g.scale((vv - 1).sqrt(), g.rot(v))))), z
    )


def prep(sa, sb):
    sm = degree((sa + sb) / 2)
    sr = degree((sb - sa) / 2)
    sall = degree(sa).union(degree(sb))
    # Enclose each sigma coefficient of D by its midpoint and the next derivative.
    dc = slope(arb_series([sm, 1], 5))
    dw = slope(arb_series([sall, 1], 5))

    def narrowed(k, comp):
        return (
            dc[comp][k]
            + sr * abs((k + 1) * dw[comp][k + 1]).union(-abs((k + 1) * dw[comp][k + 1]))
        ).intersection(dw[comp][k])

    return [
        [SigmaJet([dc[j][k] for k in range(4)]) for j in range(2)],
        [SigmaJet([narrowed(k, j) for k in range(4)]) for j in range(2)],
        sr,
    ]


def evalq(D, ww):
    """The first four sigma coefficients of q, each retaining seven w-coefficients."""
    w = SigmaJet(arb_series([ww, 1], 7))
    x = g.plus(z, g.scale(w, D))

    def inv(normal, e):
        off = g.dot(x, normal) - h
        return g.minus(
            g.minus(x, g.scale(off, normal)), g.scale((1 - off**2).sqrt(), e)
        )

    ll, nt = inv(n[5], el), inv(n[4], en)
    t = ll[0] + h
    j = (t + rt * (1 - t * t).sqrt()) / 2 - h
    return g.circle_pair((j, (j - 1) / rt), g.minus(g.scale(2, mid), nt), -1)


def boxq(pre, w0, w1):
    """First and second w-derivatives throughout the two-dimensional parameter cell."""
    dc, ds, sr = pre
    wm, wr = g.exact((w0 + w1) / 2), g.exact((w1 - w0) / 2)
    wb = g.exact(w0).union(g.exact(w1))
    q0 = evalq(dc, wm)
    qw = evalq(dc, wb)
    qs = evalq(ds, wb)
    out = []
    for k in (1, 2):
        vs = []
        for j in range(2):
            # Outer Taylor in sigma at fixed w, inner Taylor in w for the retained coefficients.
            r = sr**3 * abs(qs[j][3][k])
            for i in range(3):
                er = comb(k + 4, 4) * wr**4 * abs(qw[j][i][k + 4])
                er += sum(
                    comb(k + l, l) * wr**l * abs(q0[j][i][k + l]) for l in range(1, 4)
                )
                if i:
                    er += abs(q0[j][i][k])
                r += sr**i * er
            vs.append(factorial(k) * (q0[j][0][k] + (-r).union(r)))
        out.append(vs)
    return out


def verify():
    ctx.prec = 160
    ctx.cap = 10
    count_sigma = 100
    count_w = 40
    low_sigma, high_sigma = F("1.50"), F("1.60")
    low_w, high_w = F(".8"), F(1)
    min_curvature_left = []
    max_curvature_right = []
    sharp_lows = []
    sharp_highs = []
    lower_losses = []
    upper_losses = []
    graph_lows = []
    graph_highs = []
    speed_lows = []
    speed_highs = []
    for j in range(count_sigma):
        sa, sb = (
            low_sigma + (high_sigma - low_sigma) * i / count_sigma for i in (j, j + 1)
        )
        bounds = prep(sa, sb)
        for param, target, endpoint in (
            (low_w, min_curvature_left, "left"),
            (high_w, max_curvature_right, "right"),
        ):
            vel, acc = boxq(bounds, param, param)
            kappa = g.cross(vel, acc) / g.norm(vel) ** 3
            target.append(kappa.lower() if endpoint == "left" else kappa.upper())
        low = []
        high = []
        gl = []
        gh = []
        density_low = []
        density_high = []
        for k in range(count_w):
            a, b = (low_w + (high_w - low_w) * i / count_w for i in (k, k + 1))
            vel, acc = boxq(bounds, a, b)
            speed = -vel[0]
            det = g.cross(vel, acc)
            assert speed > 0 and det > 0
            density = det.root(3)
            density_low.append(density.lower())
            density_high.append(density.upper())
            concavity = det / speed**3
            low.append(speed.lower())
            high.append(speed.upper())
            gl.append(concavity.lower())
            gh.append(concavity.upper())
        width = g.exact((high_w - low_w) / count_w)
        sharp_lows.append(((sum(density_low) * width) ** 3 / 12).lower())
        sharp_highs.append(((sum(density_high) * width) ** 3 / 12).upper())
        length_low = sum(low) * width
        length_high = sum(high) * width
        lower_losses.append((min(gl) * length_low**3 / 12).lower())
        upper_losses.append((max(gh) * length_high**3 / 12).upper())
        graph_lows.append(min(gl))
        graph_highs.append(max(gh))
        speed_lows.append(min(low))
        speed_highs.append(max(high))
    assert min(min_curvature_left) > max(max_curvature_right) > 0
    return {
        "sigma_degrees": ["1.50", "1.60"],
        "paired_parameters": [".8", "1"],
        "sigma_cells": count_sigma,
        "parameter_cells_per_sigma": count_w,
        "arb_precision_bits": ctx.prec,
        "left_curvature_lower": str(min(min_curvature_left)),
        "right_curvature_upper": str(max(max_curvature_right)),
        "endpoint_curvature_separation_lower": str(
            (min(min_curvature_left) - max(max_curvature_right)).lower()
        ),
        "graph_concavity_lower": str(min(graph_lows)),
        "graph_concavity_upper": str(max(graph_highs)),
        "horizontal_speed_lower": str(min(speed_lows)),
        "horizontal_speed_upper": str(max(speed_highs)),
        "asymptotic_sharp_coefficient_lower": str(min(sharp_lows)),
        "asymptotic_sharp_coefficient_upper": str(max(sharp_highs)),
        "N_squared_minimum_area_loss_lower": str(min(lower_losses)),
        "N_squared_equal_horizontal_mesh_area_loss_upper": str(max(upper_losses)),
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("output", type=Path)
    args = cli.parse_args()
    args.output.write_text(json.dumps(verify(), indent=2) + "\n")
