"""Transfer the R221 full tree using per-leaf curve walls and three simultaneous point budgets."""

import argparse
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from math import isqrt
from pathlib import Path
import importlib.util
import json
import sys
import tarfile

import flint
from flint import arb, ctx
import numpy as np

S = F(".999999006")
EPS = 1 - S
R = F(".703")
r = F(1, 2)
sr = S * r
GOAL = F(120689264254125501, 2**57)
M5 = F(".8372950062210252434563")
A0 = F(".4742688878")
M0 = F(57663128211, 124250000000)
KNOTS = [(A0, M0), (F(".49"), F(".48509")), (F(".495"), F(".49119")), (r, F(".49628"))]
Q = 10**8


def ar(q):
    z = F(q)
    return arb(z.numerator) / z.denominator


def floor_scaled(v):
    return int((v.lower() * Q).floor().unique_fmpz())


def sqrt_up(q):
    den = 10**24
    lo = isqrt(q.numerator * den**2 // q.denominator)
    return F(lo + int(lo * lo * q.denominator < q.numerator * den**2), den)


def cap_knot(a, m):
    assert 0 < m < a <= r
    v = (1 + S) * r * r + EPS * m * m
    budget = 4 * S * S * (a - m)
    upper = R + EPS * max(z for _, z in KNOTS)
    # Cover a common (x,e) domain for all rates, so convexity in m interpolates exact knots.
    leaves = []

    def check(lo, hi, depth):
        assert depth < 70
        numerator = v - 2 * m * lo
        if numerator >= 0:
            den = lo + r
        else:
            den = hi + r
        # e<=x-r when the multiplier of e is positive; otherwise its product is <=0.
        bound = max(F(0), (1 + S) * (hi - r) + numerator / den)
        if lo > r:
            den = (lo * lo - r * r) if numerator >= 0 else (hi * hi - r * r)
            H = sqrt_up((GOAL - S * S * M5) * sqrt_up(hi * hi - sr * sr) / sr)
            bound = min(bound, H * max(F(0), 1 + S + numerator / den))
        if bound < budget:
            leaves.append((lo, hi, budget - bound))
            return
        mid = (lo + hi) / 2
        check(lo, mid, depth + 1)
        check(mid, hi, depth + 1)

    check(r, upper, 0)
    assert leaves[0][0] == r and leaves[-1][1] == upper
    assert all(u[1] == v[0] for u, v in zip(leaves, leaves[1:]))
    return dict(
        wall=str(a),
        rate=str(m),
        slabs=len(leaves),
        slack=str(min(z for _, _, z in leaves)),
    )


def read_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def main(comp, tree, out, rho):
    assert __debug__
    ctx.prec = 160
    pi = arb.pi()
    sys.path.insert(0, str(comp))
    from certify_tomar_six_high import SourceTree
    from certify_finite_point_leaf import point_tables

    br = read_module(
        "signed_row_peak", comp / "balanced-cap-transfer/certify_balanced_rows.py"
    )
    hashes = {}

    def read(name):
        raw = (comp / name).read_bytes()
        hashes[name] = sha256(raw).hexdigest()
        return json.loads(raw)

    old = read("balanced-cap-transfer/verified.json")
    base = read("six-input-083745/points/normal-gap-083745-audit.json")
    original = read("projected-cap-transfer/verified.json")
    cloud = read("finite-point-leaf-input.json")
    bank = read("finite-point-global-input.json")
    assert F(base["wall"]) == A0 and F(old["reference_projection_rate"]) == M0
    assert F(original["five_minimum_lower"]) == M5 and F(original["radius_upper"]) == F(
        ".528617"
    )
    assert all(F(original[k]) == z for k, z in [("scale", S), ("threshold", GOAL)])
    # Bind the generated guard, source tree decoder, and angle table to R233.
    source = old["angle_partition"]["source"]
    for name, expected in source["hashes"].items():
        hashes[name] = sha256((comp / name).read_bytes()).hexdigest()
        assert hashes[name] == expected
    for name in [
        "balanced-cap-transfer/certify_balanced_rows.py",
        "projected-cap-transfer/verify.py",
    ]:
        hashes[name] = sha256((comp / name).read_bytes()).hexdigest()
    assert (
        hashes["balanced-cap-transfer/certify_balanced_rows.py"]
        == old["source_sha256"]["certify_balanced_rows.py"]
    )
    assert (
        hashes["projected-cap-transfer/verify.py"] == old["source_sha256"]["verify.py"]
    )
    assert (
        hashes["finite-point-leaf-input.json"]
        == old["angle_partition"]["source"]["hashes"]["finite-point-leaf-input.json"]
    )
    assert (
        hashes["finite-point-global-input.json"]
        == old["angle_partition"]["source"]["hashes"]["finite-point-global-input.json"]
    )
    archive = comp / "six-input-083745/volumes/native45-receipts.tar.gz"
    assert sha256(archive.read_bytes()).hexdigest() == source["archive_sha256"]
    with tarfile.open(archive) as z:
        raw_scope = z.extractfile("inputs/ordinary/scope.json").read()
    assert (
        sha256(raw_scope).hexdigest()
        == source["archive_members_sha256"]["inputs/ordinary/scope.json"]
    )
    scope = json.loads(raw_scope)
    keys = scope["leaves"] + scope["deferred"]
    assert len(keys) == len(set(map(tuple, keys))) == 3273
    tr = SourceTree(tree)
    boxes = [tr.leaf(*key) for key in keys]
    ids = sorted({i for cell in bank["cells"] for row in cell for i in row[:3]})
    V = cloud["directions"]
    norms = [arb(x * x + y * y).sqrt() for x, y in V]
    U = [(x / l, y / l) for (x, y), l in zip(V, norms)]
    angles = [arb.atan2(y, x) for x, y in V]

    @lru_cache(None)
    def support(n, i, span):
        half = pi / n
        center = 1 / (2 * (half / 2).cos())
        mid, rad = ar((span[0] + span[1]) / 2), ar((span[1] - span[0]) / 2)
        vals = [abs(angles[i] - mid - 2 * j * half) for j in range(-n, n + 1)]
        d0 = min(x.lower() for x in vals).union(min(x.upper() for x in vals))
        dl = max(arb(0), (d0 - rad).lower())
        dh = min(half.upper(), (d0 + rad).upper())

        def H(d):
            cand = []
            if not d > half / 2:
                cand.append(center * d.cos())
            if not d < half / 2:
                cand.append(1 - center * (half - d).cos())
            assert cand
            return cand[0] if len(cand) == 1 else cand[0].union(cand[1])

        return H(dh).lower(), H(dl).upper()

    @lru_cache(None)
    def shift(i, k, span):
        u = U[i][k]
        l, h = span
        return max((u * ar(l)).upper(), (u * ar(h)).upper())

    @lru_cache(None)
    def wall(n, i, angle, tx, ty):
        hl, hh = support(n, i, angle)
        T = shift(i, 0, tx) + shift(i, 1, ty)
        minwall = max(
            hl, (ar(r) - T).lower(), ((ar(r) - T - ar(S) * hh) / ar(EPS)).lower()
        )
        return min(Q // 2, max(int(A0 * Q), floor_scaled(minwall)))

    @lru_cache(None)
    def rate(ticks):
        x = max(A0, F(int(ticks), Q))
        for (a, m), (b, p) in zip(KNOTS, KNOTS[1:]):
            if x <= b:
                return int(Q * (m + (p - m) * (x - a) / (b - a)))
        raise AssertionError

    rates = np.zeros((len(boxes), len(V)), dtype=np.int64)
    for j, b in enumerate(boxes):
        for i in ids:
            ticks = min(
                wall(
                    n,
                    i,
                    tuple(b[angle]) if angle is not None else (F(0), F(0)),
                    tuple(b[t]),
                    tuple(b[t + 1]),
                )
                for n, angle, t in ((3, None, 0), (5, 2, 3), (7, 5, 6), (9, 8, 9))
            )
            rates[j, i] = rate(ticks)
    assert np.max(rates) <= Q // 2 and np.min(rates) >= 0 and 3 * Q * Q < 2**63
    _, ptr, flat = point_tables({**cloud, "cells": bank["cells"]})
    P = cloud["points"]
    dP = cloud["denominator"]
    arc = [((2 * pi * j / 384).cos(), (2 * pi * j / 384).sin()) for j in range(385)]

    @lru_cache(None)
    def joint(ids, weights):
        den = sum(w * norms[i] for i, w in zip(ids, weights))
        coef = np.array(
            [floor_scaled(w * norms[i] / den) for i, w in zip(ids, weights)],
            dtype=np.int64,
        )
        assert all(0 <= x <= Q for x in coef)
        values = rates[:, ids] @ coef
        j = int(np.argmin(values))
        return F(int(values[j]), Q * Q), j

    worst = None
    intersections = 0
    for c, items in enumerate(bank["cells"]):
        for j, row in enumerate(items):
            inds, points = tuple(row[:3]), row[3:]
            packed = flat[8 * (ptr[c] + j) : 8 * (ptr[c] + j + 1)]
            assert list(inds) == packed[:3]
            wt = tuple(packed[3:6])
            assert all(w > 0 for w in wt) and all(
                sum(w * V[i][k] for i, w in zip(inds, wt)) == 0 for k in (0, 1)
            )
            den = sum(w * norms[i] for i, w in zip(inds, wt)) * dP
            A = sum(
                w * (V[i][0] * P[p][0] + V[i][1] * P[p][1])
                for w, i, p in zip(wt, inds, points)
            )
            B = sum(
                w * (V[i][1] * P[p][0] - V[i][0] * P[p][1])
                for w, i, p in zip(wt, inds, points)
            )
            _, up, _ = br.peak(A, B, arc[c], arc[c + 1])
            q = (up / den).upper()
            assert q > 0
            budget, idx = joint(inds, wt)
            safe = (ar(EPS * budget) / q).lower()
            assert safe > ar(rho), (c, j, keys[idx], safe, rho)
            if worst is None or safe < worst[0]:
                worst = (safe, c, j, keys[idx], str(budget), str(q))
            intersections += len(keys)
    expected = len(keys) * sum(map(len, bank["cells"]))
    assert intersections == expected == 3273 * 2493
    cap = [cap_knot(a, m) for a, m in KNOTS]
    dp = F(original["point_diameter_squared"])
    factor = 1 / (dp * (1 - rho) ** 2)
    lower = factor * GOAL
    assert lower > F(old["transfer_area_lower"]) and S * S * factor < 1
    pl = 24
    val = int(lower * 10**pl)
    result = dict(
        complete=True,
        scale=str(S),
        threshold=str(GOAL),
        point_diameter_squared=str(dp),
        point_shrink=str(rho),
        lower=str(lower),
        factor_squared=str(factor),
        strict_floor=f"{val//10**pl}.{val%10**pl:0{pl}d}",
        improvement=str(lower - F(old["transfer_area_lower"])),
        minimum_allowed_shrink=str(worst[0]),
        attaining_angle_cell=worst[1],
        attaining_row=worst[2],
        attaining_leaf=worst[3],
        weighted_rate_lower=worst[4],
        selected_projection_upper=worst[5],
        joint_leaf_row_pairs=intersections,
        distinct_triples=joint.cache_info().currsize,
        cap_knots=cap,
        native_original_leaves=len(keys),
        quantization_denominator=Q,
        tree_manifest_sha256=sha256((tree / "SHA256SUMS").read_bytes()).hexdigest(),
        verifier_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        archive_sha256=sha256(archive.read_bytes()).hexdigest(),
        source_hashes=hashes,
        python_flint=flint.__version__,
        precision_bits=ctx.prec,
    )
    out.write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")
    print(result["strict_floor"], str(worst), cap, flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("computation", type=Path)
    p.add_argument("tree", type=Path)
    p.add_argument("out", type=Path)
    p.add_argument("--rho", type=F, default=F(".00000089591"))
    a = p.parse_args()
    main(a.computation.resolve(), a.tree.resolve(), a.out, a.rho)
