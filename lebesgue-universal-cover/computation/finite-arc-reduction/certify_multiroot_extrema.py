"""Reconstruct polarized Reuleaux polygons for a fixed two-cap host with many genuine roots."""

import argparse
from bisect import bisect_right
from fractions import Fraction as Q
from hashlib import sha256
import json
import math
from pathlib import Path
import platform

import flint
from flint import arb, ctx


def number(q):
    return arb(q.numerator) / q.denominator


def unit(t):
    sine, cosine = t.sin_cos()
    return cosine, sine


def perp(v):
    return -v[1], v[0]


def add(u, v):
    return tuple(x + y for x, y in zip(u, v))


def sub(u, v):
    return tuple(x - y for x, y in zip(u, v))


def scale(v, c):
    return tuple(x * c for x in v)


def dot(u, v):
    return sum(x * y for x, y in zip(u, v))


def reconstruct(r):
    assert r % 2 == 1 and math.gcd(3 * r, 5) == 1
    n = 3 * r
    a = arb(1) / (8 * (n * n - 1))
    b = a / 100
    pi = arb.pi()
    L = Q(1, 3)
    alpha = Q(1, 18)
    angle = lambda q: pi * number(q)
    roots = [Q(j, n) for j in range(r)]

    def h(t):
        return arb(1) / 2 + a * (n * t).sin() + b * ((5 * t).sin() + t.sin())

    def point(t):
        u = unit(t)
        dh = n * a * (n * t).cos() + b * (5 * (5 * t).cos() + t.cos())
        return add(scale(u, h(t)), scale(perp(u), dh))

    mesh = {Q(0), Q(1)}
    for root in roots:
        for j in range(3):
            for offset in (Q(0), -alpha, alpha):
                mesh.add((root + j * L + offset) % 1)
    mesh = sorted(mesh)
    nodes = {q: h(angle(q)) for q in mesh}
    pieces = []
    prev = None
    switches = 0
    first = None
    for l, rr in zip(mesh, mesh[1:]):
        mid = (l + rr) / 2
        sector = int(mid / L)
        sign = int(((mid % L) * n)) % 2 == 0
        upper = sign == (sector % 2 == 0)
        left, right = angle(l), angle(rr)
        pl, pr = point(left), point(right)
        chord = sub(pr, pl)
        if upper:
            chord = sub(sub(unit(right), unit(left)), chord)
        center_angle = arb.atan2(-chord[0], chord[1])
        ell = (dot(chord, chord).sqrt() / 2).asin()
        cuts = (left, center_angle - ell, center_angle + ell, right)
        assert all(x < y for x, y in zip(cuts, cuts[1:]))
        local = []
        p = pl
        for j, (start, end) in enumerate(zip(cuts, cuts[1:])):
            radius = int(upper) if j != 1 else int(not upper)
            if first is None:
                first = radius
            if prev is not None and prev != radius:
                switches += 1
            prev = radius
            c = sub(p, scale(unit(start), radius))
            p = add(c, scale(unit(end), radius))
            local.append((start, end, radius, c))
        assert all(s.contains(0) for s in sub(p, pr))
        pieces.append(local)
    arcs = switches + int(prev != 1 - first)
    assert 0 < arcs <= 33 * r

    def fitted(t):
        if t in nodes:
            return nodes[t]
        ind = bisect_right(mesh, t) - 1
        v = angle(t)
        fits = [
            (radius, c) for left, right, radius, c in pieces[ind] if left < v < right
        ]
        assert len(fits) == 1, (t, ind)
        radius, c = fits[0]
        return radius + dot(c, unit(v))

    # Roots and pose translations depend only on preserved hexagon records.
    rows = []
    for root in roots:
        dx = arb(1) / 2 - nodes[root]
        dy = (arb(1) / 2 - nodes[root + L] - dx / 2) / (arb(3).sqrt() / 2)
        assert (
            nodes[root] - nodes[root + L] + nodes[root + 2 * L] - arb(1) / 2
        ).contains(0)
        for k in range(6):
            for sign in (1, -1):
                gaps = []
                for cap in (alpha, -alpha):
                    f = (sign * cap + k * L) % 2
                    src = (root + f) % 2
                    value = nodes[src] if src < 1 else 1 - nodes[src - 1]
                    gaps.append(value + dot((dx, dy), unit(angle(f))) - arb(1) / 2)
                gap = max(g.lower() for g in gaps)
                assert gap > 0 and gap > 43 * a / 100
                rows.append(gap)
    # Check the ordering between old and new F, separately on every sign interval.
    min_gain = None
    samples = 0
    for j in range(r):
        for k in range(1, 11):
            t = Q(j, n) + Q(13 * k, 143 * n)
            original = 3 * a * (n * angle(t)).sin()
            new = fitted(t) - fitted(t + L) + fitted(t + 2 * L) - arb(1) / 2
            gain = (-1) ** j * (new - original)
            assert gain > 0 or gain.contains(0)
            assert (-1) ** j * new > 0
            min_gain = (
                gain.lower()
                if min_gain is None or gain.lower() < min_gain
                else min_gain
            )
            samples += 1
    return {
        "n": n,
        "root_classes": r,
        "distinct_hexagon_placements_lower": (r + 1) // 2,
        "half_circle_cells": len(mesh) - 1,
        "reconstructed_unit_arcs": arcs,
        "pose_tests": len(rows),
        "min_pose_violation": str(min(rows)),
        "crosscheck_angle_samples": samples,
        "sample_minimum_outward_F_change": str(min_gain),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    ctx.prec = 192
    rows = [reconstruct(r) for r in (3, 7, 11, 21, 31)]
    result = {
        "precision_bits": ctx.prec,
        "python": platform.python_version(),
        "python_flint": flint.__version__,
        "script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "host_caps_pi_units": ["-1/18", "1/18"],
        "circles": rows,
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(
        [
            (x["root_classes"], x["reconstructed_unit_arcs"], x["pose_tests"])
            for x in rows
        ]
    )
