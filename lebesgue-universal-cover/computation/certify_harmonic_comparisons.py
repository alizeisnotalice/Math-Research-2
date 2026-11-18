"""Certify a cubic inside the Round 56 hull and a common upper for eight moving bodies."""

from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
from flint import arb, ctx
from certify_nine import (
    PI_QUALITY,
    ar,
    append_point,
    arc_body,
    fourier_body,
    FOURIER_TERMS,
    hull,
    arithmetic,
)
from harmonic_inner_material import polygon, DEN

ctx.prec = PI_QUALITY
HERE = Path(__file__).parent


def cubic(source):
    rho, tx, ty = map(ar, [".00310264", ".0006235264", "-.0089536526"])
    sn, cs = (3 * rho).sin_cos()
    margin = Q(".0000002")
    p = polygon(source, 256)
    for j, (x, y) in enumerate(p):
        a, b = p[(j + 1) % len(p)]
        X, Y = b - y, x - a
        L = arb(X * X + Y * Y).sqrt()
        c, s = arb(X) / L, arb(Y) / L
        support = (
            ar(Q(1, 2))
            + ((4 * c**3 - 3 * c) * cs + (3 * s - 4 * s**3) * sn) / 16
            + tx * c
            + ty * s
        )
        assert arb(X * x + Y * y) / (L * DEN) - support > ar(margin)
    return dict(
        angle=".00310264",
        translation=[".0006235264", "-.0089536526"],
        interior_disk_padding=str(margin),
        degree=256,
        strict_edge_count=len(p),
    )


def upper(params, k):
    x = [ar(z) for z in params]
    assert len(x) == 20
    pi = arb.pi()
    pts = []
    for j in range(8 * k):
        s, c = (2 * j * pi / (8 * k)).sin_cos()
        append_point(pts, c / 2, s / 2)
    arc_body(pts, k, pi, 3, arb(0), *x[:2])
    for n, j in ((5, 2), (7, 5), (9, 8)):
        arc_body(pts, k, pi, n, *x[j : j + 3])
    for j, terms in enumerate(FOURIER_TERMS + [((11, Q(1, 240)),)]):
        i = 11 + 3 * j
        fourier_body(pts, (8 if j == 0 else 16) * k, pi, terms, *x[i : i + 3])
    # R=h+h'' is between zero and one, |R'|<=11/2 for the added body.
    assert 1 + Q(11, 2) ** 2 < Q(32, 5) ** 2
    return arithmetic(hull(pts), k)


if __name__ == "__main__":
    p = (HERE / "056-free-seven.json").read_bytes()
    r = (HERE / "harmonic-eight-placement.json").read_bytes()
    out = upper(json.loads(r)["best"]["params"], 1024)
    fixed = Q(json.loads((HERE / "harmonic-strong-lower.json").read_bytes())["target"])
    seven = Q(json.loads((HERE / "056-free-area.json").read_bytes())["area"]["upper"])
    free = Q(out["upper"])
    assert fixed == Q(".83828") and free < Q(".838250180")
    coord = Q(".0000055")
    delta = coord * Q(19, 12)
    assert Q(2) < (Q(19, 12) - Q(1, 6)) ** 2
    assert ar(Q(1, 2)) * (arb.pi() / 10).tan() < ar(Q(1, 6))
    assert all(sum(n * abs(a) for n, a in terms) < Q(1, 6) for terms in FOURIER_TERMS)
    body_dist = Q(".0000044")
    assert fixed - seven > Q(".00031768") and fixed - free > Q(".00002982")
    assert fixed / (1 + 2 * delta) ** 2 - free > Q(".00000062")
    assert fixed / (1 + 2 * body_dist) ** 2 - free * (1 + 2 * body_dist) ** 2 > Q(
        ".00000031"
    )
    comparisons = dict(
        fixed_minimum_strict_lower=str(fixed),
        free_upper=str(free),
        fixed_addition_gap_lower=str(fixed - seven),
        fixed_minus_free_lower=str(fixed - free),
        closed_base_parameter_halfwidth=str(coord),
        robust_fixed_minus_free_lower=str(fixed / (1 + 2 * delta) ** 2 - free),
        tested_body_hausdorff_radius=str(body_dist),
        perturbed_body_gap_lower=str(
            fixed / (1 + 2 * body_dist) ** 2 - free * (1 + 2 * body_dist) ** 2
        ),
    )
    print(
        json.dumps(
            dict(
                source_sha256=sha256(p).hexdigest(),
                joint_sha256=sha256(r).hexdigest(),
                cubic=cubic(p),
                joint_degree=1024,
                eight=out,
                comparisons=comparisons,
            ),
            indent=2,
        )
    )
