"""Check the exact replacement P30, all angular increases, and native integer score slots."""

from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import json, sys
from flint import arb, ctx


def load(path):
    return json.loads(path.read_bytes())


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def ar(x):
    f = F(x)
    return arb(f.numerator) / f.denominator


def run(comp, f, out):
    assert __debug__
    ctx.prec = 192
    sys.path.insert(0, str(comp))
    from certify_finite_point_leaf import point_tables, Q, U

    pins = {
        "finite-point-leaf-input.json": "b99b3bf66adf904b8b99ba1e1c7f8e75abc1ea6d8e00728b9553c20b961ffebc",
        "finite-point-global-input.json": "348f983e915e77a50b55d627d3d2728a99377f7283dd67848896b0bcfff5469f",
        "certify_finite_point_leaf.py": "f2992ef5ec1ab1675a25a3d4e8321cd7b20c940c9da7a918821eb8e0e1a8d2ab",
    }
    assert all(digest(comp / k) == h for k, h in pins.items())
    assert (
        digest(f) == "1afef7eabda436e4950853fc4623aa2b3956b625b747bdd7053732fecfe5eb58"
    )
    P = load(comp / "finite-point-leaf-input.json")
    Y = load(f)
    bank = load(comp / "finite-point-global-input.json")["cells"]
    assert (
        P["direction_scale"] == Y["direction_scale"] == U
        and P["directions"] == Y["directions"]
    )
    V = P["directions"]
    dp = P["denominator"]
    dy = Y["denominator"]
    Z = [
        (x * dp - a * dy, y * dp - b * dy)
        for (x, y), (a, b) in zip(Y["points"], P["points"], strict=True)
    ]
    pdiam, ptr, flat = point_tables({**P, "cells": bank})
    ydiam, ptr1, new = point_tables({**Y, "cells": bank})
    assert ptr == ptr1 and len(Z) == 30 and ydiam < dy * dy
    for k in range(7):
        assert flat[k::8] == new[k::8]
    deltas = [b - a for a, b in zip(flat[7::8], new[7::8])]
    angles = [(2 * arb.pi() * i / 384).sin_cos() for i in range(385)]
    min_actual = None
    for c, cell in enumerate(bank):
        for j, row in enumerate(cell):
            k = ptr[c] + j
            ids, pts = row[:3], row[3:]
            w = flat[8 * k + 3 : 8 * k + 6]
            assert ids == flat[8 * k : 8 * k + 3] and sum(w) == flat[8 * k + 6]
            assert all(sum(wi * V[n][d] for wi, n in zip(w, ids)) == 0 for d in (0, 1))
            den = dp * dy * U * sum(w)
            A = sum(
                h * (V[n][0] * Z[p][0] + V[n][1] * Z[p][1])
                for h, n, p in zip(w, ids, pts)
            )
            B = sum(
                h * (V[n][1] * Z[p][0] - V[n][0] * Z[p][1])
                for h, n, p in zip(w, ids, pts)
            )
            # The positive endpoint values imply positivity and concavity on a cell of width <pi.
            v = min(((A * cs + B * sn) / den).lower() for sn, cs in angles[c : c + 2])
            assert v > ar(".0002549")
            if min_actual is None or v < min_actual[0]:
                min_actual = (v, c, j)
    assert (
        len(bank) == 384
        and len(deltas) == 2493
        and min(deltas) > int(F(".0002549") * Q)
    )
    arg = deltas.index(min(deltas))
    c = next(j for j, (u, v) in enumerate(zip(ptr, ptr[1:])) if u <= arg < v)
    result = dict(
        complete=True,
        inputs_sha256=pins,
        replacement_sha256=digest(f),
        verifier_sha256=digest(Path(__file__)),
        diameter_squared=str(F(ydiam, dy**2)),
        integer_diameter_squared=ydiam,
        diameter_pairs=435,
        normalized_projection_gain_lower=str(min_actual[0]),
        attaining_angle_cell=min_actual[1],
        attaining_row=min_actual[2],
        all_original_material_indices_preserved=True,
        directions=len(V),
        cells=len(bank),
        rows=len(deltas),
        distinct_triples=len({tuple(row[:3]) for cell in bank for row in cell}),
        native_score_scale=Q,
        minimum_native_score_gain=min(deltas),
        maximum_native_score_gain=max(deltas),
        minimum_native_gain_cell=c,
        minimum_native_gain_row=arg - ptr[c],
        native_normalization_spatial=U,
        reflected=True,
    )
    out.write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")
    print(
        result["normalized_projection_gain_lower"],
        min(deltas),
        max(deltas),
        min_actual[1:],
        flush=True,
    )


if __name__ == "__main__":
    a = ArgumentParser(description=__doc__)
    a.add_argument("computation", type=Path)
    a.add_argument("replacement", type=Path)
    a.add_argument("out", type=Path)
    p = a.parse_args()
    run(p.computation.resolve(), p.replacement, p.out)
