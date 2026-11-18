"""Certify a non-similar deformation against all closed R233 angular rows."""

from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import lcm
from pathlib import Path
import importlib.util
import json
import sys
import flint
from flint import arb, ctx

# Only the point numbers in the upper closed half-plane are needed.
# Entries are ten times the displacement divided by epsilon.
POSITIVES = {
    0: (-11, 0),
    2: (-15, -25),
    4: (-25, -21),
    6: (-13, 25),
    8: (-25, -14),
    10: (25, 17),
    12: (-25, -6),
    14: (-25, -19),
    16: (25, 25),
    18: (25, 25),
    20: (25, 25),
    22: (-21, 25),
    24: (-25, 25),
    26: (8, -25),
    28: (-4, 25),
    29: (0, 0),
}
S = F(".999999006")
EPS = 1 - S
G = F(120689264254125501, 2**57)
M0 = F(57663128211, 124250000000)


def ar(q):
    q = F(q)
    return arb(q.numerator) / q.denominator


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def dump(path, obj):
    path.write_text(json.dumps(obj, sort_keys=True, indent=2) + "\n")


def run(comp, dest):
    assert __debug__
    ctx.prec = 192
    br = comp / "balanced-cap-transfer" / "certify_balanced_rows.py"
    oldpath = comp / "balanced-cap-transfer" / "verified.json"
    assert (
        sha(oldpath)
        == "6407c817e042fb5b2649d82a2114f086b02531e7e097ef13911b8774c13ccded"
    )
    old = json.loads(oldpath.read_bytes())
    assert sha(br) == old["source_sha256"]["certify_balanced_rows.py"]
    spec = importlib.util.spec_from_file_location("balanced", br)
    b = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(b)
    hashes = {
        str(br.relative_to(comp)): sha(br),
        str(oldpath.relative_to(comp)): sha(oldpath),
    }
    for name in [
        "finite-point-leaf-input.json",
        "finite-point-global-input.json",
        "certify_finite_point_leaf.py",
    ]:
        hashes[name] = sha(comp / name)
        assert hashes[name] == old["angle_partition"]["source"]["hashes"][name]
    pp = comp / "projected-cap-transfer" / "verified.json"
    origin = json.loads(pp.read_bytes())
    assert sha(pp) == old["source_sha256"]["verified.json"]
    assert F(old["reference_projection_rate"]) == F(origin["projection_rate"]) == M0
    assert F(origin["scale"]) == S and F(origin["threshold"]) == G
    assert EPS > 0 and F(".419") < M0
    cloud = json.loads((comp / "finite-point-leaf-input.json").read_bytes())
    cells = json.loads((comp / "finite-point-global-input.json").read_bytes())["cells"]
    P = cloud["points"]
    V = cloud["directions"]
    D = cloud["denominator"]
    Z = []
    for j, (x, y) in enumerate(P):
        i = P.index([x, abs(y)])
        zx, zy = POSITIVES[i]
        Z.append((zx, zy if y >= 0 else -zy))
    assert set(POSITIVES) == {j for j, (x, y) in enumerate(P) if y >= 0}
    assert all((a, 0) == (a, b) for (_, y), (a, b) in zip(P, Z) if y == 0)
    E = lcm(D, (EPS / 10).denominator)
    Y = [
        tuple(p * int(E // D) - int(E * EPS / 10) * z for p, z in zip(t, v))
        for t, v in zip(P, Z)
    ]

    def diam(points):
        return max(
            ((x - u) ** 2 + (y - v) ** 2, i, j)
            for (i, (x, y)), (j, (u, v)) in combinations(enumerate(points), 2)
        )

    sq, i, j = diam(Y)
    d2 = F(sq, E * E)
    assert d2 < S * S
    assert {(x, -y) for x, y in Y} == set(Y)
    # Exact diameter-one similarity, dividing by the complex diametral vector.
    a, beta = Y[i][0] - Y[j][0], Y[i][1] - Y[j][1]
    norm = [(a * x + beta * y, a * y - beta * x) for x, y in Y]
    assert diam(norm)[0] == sq * sq
    sys.path.insert(0, str(comp))
    from certify_finite_point_leaf import point_tables

    _, ptr, flat = point_tables({**cloud, "cells": cells})
    arc = [
        ((2 * arb.pi() * j / 384).cos(), (2 * arb.pi() * j / 384).sin())
        for j in range(385)
    ]
    lengths = [arb(x * x + y * y).sqrt() for x, y in V]
    peak = None
    for cell, items in enumerate(cells):
        for pos, row in enumerate(items):
            inds, ps = row[:3], row[3:]
            f = flat[8 * (ptr[cell] + pos) : 8 * (ptr[cell] + pos + 1)]
            w = f[3:6]
            assert inds == f[:3] and all(z > 0 for z in w)
            assert all(sum(h * V[k][j] for h, k in zip(w, inds)) == 0 for j in range(2))
            den = 10 * sum(h * lengths[k] for h, k in zip(w, inds))
            A = sum(
                h * (V[k][0] * Z[t][0] + V[k][1] * Z[t][1])
                for h, k, t in zip(w, inds, ps)
            )
            B = sum(
                h * (V[k][1] * Z[t][0] - V[k][0] * Z[t][1])
                for h, k, t in zip(w, inds, ps)
            )
            lo, hi, kind = b.peak(A, B, arc[cell], arc[cell + 1])
            assert hi / den < ar(".419"), (cell, pos, hi / den)
            if peak is None or (hi / den).upper() > peak[0]:
                peak = ((hi / den).upper(), (lo / den).lower(), cell, pos, kind)
    assert sum(map(len, cells)) == 2493 and len(cells) == 384
    dest.mkdir(parents=True, exist_ok=True)
    point_file = dest / "points.json"
    dump(
        point_file,
        dict(
            original_point_order=True,
            deformation_numerators=Z,
            deformation_denominator=10,
            deformed=dict(points=Y, denominator=E),
            exact_unit=dict(points=norm, denominator=sq),
        ),
    )
    low = G / (S * S)
    floor = (low.numerator * 10**24) // low.denominator
    result = dict(
        complete=True,
        source_sha256=hashes,
        verifier_sha256=sha(Path(__file__)),
        points_sha256=sha(point_file),
        scale=str(S),
        threshold=str(G),
        old_rate=str(M0),
        area_lower=str(low),
        strict_floor=f"{floor//10**24}.{floor%10**24:024d}",
        angular_rows=2493,
        total_angle_cells=384,
        handedness="reflection verified",
        weighted_max_upper=str(peak[0]),
        weighted_max_lower=str(peak[1]),
        attaining_cell=peak[2],
        attaining_row=peak[3],
        attaining_type=peak[4],
        squared_diameter=str(d2),
        diametral_pair=[i, j],
        squared_diameter_slack=str(S * S - d2),
        all_pairs=435,
        normalized_exact_squared_diameter=1,
        precision_bits=ctx.prec,
        flint=flint.__version__,
    )
    dump(dest / "proof.json", result)
    print(result["strict_floor"], peak, i, j, float((S * S - d2)), flush=True)


if __name__ == "__main__":
    p = ArgumentParser(description=__doc__)
    p.add_argument("computation", type=Path)
    p.add_argument("destination", type=Path)
    a = p.parse_args()
    run(a.computation.resolve(), a.destination)
