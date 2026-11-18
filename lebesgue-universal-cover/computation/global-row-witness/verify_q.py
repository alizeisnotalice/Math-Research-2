"""Bracket the loss required to replace the literal Q11 support rows at diameter s."""

from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import lcm
from pathlib import Path
import json, sys
from flint import arb, ctx


def ar(x):
    f = F(x)
    return arb(f.numerator) / f.denominator


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def load(path):
    return json.loads(path.read_bytes())


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def output(path, data):
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")


def run(comp, cert, out):
    assert __debug__
    ctx.prec = 192
    sys.path.insert(0, str(comp))
    from build_q11_leaf_bank import EPS, reflected_average
    from certify_finite_point_leaf import point_tables, Q, U
    from replay_q11_with_p30 import PointProof

    S = F(".999999006")
    e = 1 - S
    P, den = reflected_average(load(comp / "multi-host-points.json"))
    bank = load(comp / "q11-leaf-bank.json")
    cells = bank["cells"]
    V = bank["directions"]
    datum = dict(**bank, points=P, denominator=den, direction_scale=U)
    native = PointProof(datum, EPS)
    dia, ptr, raw = point_tables(datum)
    flat = list(map(int, native.triples))
    assert raw[7::8] == [v + 1 for v in flat[7::8]]
    assert all(flat[j::8] == raw[j::8] for j in range(7))
    assert (
        [int(v) for v in native.pointers] == ptr
        and len(cells) == 384
        and ptr[-1] == 2258
    )
    N = len(P)
    assert (
        N == 11
        and P[0] == [0, 0]
        and [(x, -y) for x, y in P[1:]] == list(map(tuple, reversed(P[1:])))
    )
    unit = [(ar(F(x, den)), ar(F(y, den))) for x, y in P]
    tt = [(2 * arb.pi() * i / 384).sin_cos() for i in range(385)]
    rows = []
    for c, items in enumerate(cells):
        for j, row in enumerate(items):
            k = ptr[c] + j
            ids, ps = row[:3], row[3:]
            ww = flat[8 * k + 3 : 8 * k + 6]
            assert (
                ids == flat[8 * k : 8 * k + 3]
                and min(ww) > 0
                and sum(ww) == flat[8 * k + 6]
            )
            assert all(sum(h * V[v][d] for h, v in zip(ww, ids)) == 0 for d in (0, 1))
            score = F(flat[8 * k + 7], Q)
            for end in [0, 1]:
                a = [[arb(0), arb(0)] for _ in P]
                sn, cs = tt[c + end]
                for h, v, p in zip(ww, ids, ps):
                    x, y = V[v]
                    q = ar(F(h, sum(ww) * U))
                    a[p][0] += q * (x * cs + y * sn)
                    a[p][1] += q * (-x * sn + y * cs)
                baseline = (sum(dot(t, z) for t, z in zip(a, unit)) - ar(score)) / ar(e)
                assert baseline > 0
                rows.append((c, j, end, a, baseline))
    pairs = []
    for i, j in combinations(range(N), 2):
        d = [unit[i][k] - unit[j][k] for k in range(2)]
        ell = dot(d, d).sqrt()
        if ell > ar(".9999"):
            pairs.append((i, j, [v / ell for v in d], (ar(S) - ell) / ar(e)))
        else:
            assert ell < ar(".9999")
    assert len(pairs) == 11
    obj = load(cert)
    weights = obj["dual"]
    D = obj["dual_den"]
    L = sum(w for k, w in weights if k < len(rows))
    assert abs(L - D) <= 5
    residual = [[arb(0), arb(0)] for _ in P]
    lower = arb(0)
    used = []
    for k, w in weights:
        assert w > 0
        weight = ar(F(w, L))
        if k < len(rows):
            cell, j, end, v, b = rows[k]
            lower -= weight * b
            for x, src in zip(residual, v):
                for d in range(2):
                    x[d] -= weight * src[d]
            used.append(dict(cell=cell, row=j, endpoint=end, weight=str(F(w, L))))
        else:
            i, j, n, b = pairs[k - len(rows)]
            lower -= weight * b
            for d in range(2):
                residual[i][d] += weight * n[d]
                residual[j][d] -= weight * n[d]
            used.append(dict(diameter_edge=[i, j], weight=str(F(w, L))))
    # z_0=0 by translation. Reflection ties z_{11-i}=(z_{i,x},-z_{i,y}).
    rem = sum(
        abs(residual[i][0] + residual[-i][0]) + abs(residual[i][1] - residual[-i][1])
        for i in range(1, 6)
    )
    guaranteed = (lower - rem * ar(2 / e)).lower()
    assert guaranteed > ar(".50091")
    # Rational feasible data also bound the optimal loss from above.
    vectors = [[F(x, obj["z_den"]) for x in z] for z in obj["z"]]
    vectors += [[x, -y] for x, y in reversed(vectors[1:])]
    X = [
        [F(p, den) + e * (z - F(1, 10**5) * F(p, den)) for p, z in zip(point, v)]
        for point, v in zip(P, vectors)
    ]
    dd = lcm(*(z.denominator for point in X for z in point))
    Y = [[int(x * dd), int(y * dd)] for x, y in X]
    sq = max((x - u) ** 2 + (y - v) ** 2 for (x, y), (u, v) in combinations(Y, 2))
    assert F(sq, dd * dd) < S * S
    Z = [[ar((x - F(y, den)) / e) for x, y in zip(a, b)] for a, b in zip(X, P)]
    maxloss = max(
        (-sum(dot(t, z) for t, z in zip(a, Z)) - baseline).upper()
        for _, _, _, a, baseline in rows
    )
    assert maxloss < ar(".50094")
    # In each angular cell the strictly positive support sinusoid is concave.
    assert all(ar(F(flat[8 * k + 7], Q)) > ar(e) * maxloss for k in range(ptr[-1]))
    _, ptrnew, new = point_tables({**datum, "points": Y, "denominator": dd})
    assert ptrnew == ptr
    drops = [v - u for v, u in zip(flat[7::8], new[7::8])]
    assert max(drops) <= 134
    dst = out.parent
    output(
        dst / "q-squeezed.json",
        dict(points=Y, denominator=dd, directions=V, direction_scale=U),
    )
    result = dict(
        complete=True,
        rows=ptr[-1],
        cells=len(cells),
        directions=len(V),
        native_reflection_penalty_ticks=1,
        exact_native_replacement_maximum_score_drop=max(drops),
        maximum_continuous_loss_rate_upper=str(maxloss),
        all_fixed_label_clouds_necessary_loss_rate_lower=str(guaranteed),
        dual_l1_residual=str(rem),
        dual_rows_and_diameters=used,
        scale=str(S),
        squared_diameter=str(F(sq, dd * dd)),
        loss_definition="max of (frozen integer lower score minus physical projected average) / (1-s), over 2258 rows and their closed angle cells",
        source_sha256={
            name: digest(comp / name)
            for name in [
                "multi-host-points.json",
                "q11-leaf-bank.json",
                "build_q11_leaf_bank.py",
                "replay_q11_with_p30.py",
                "certify_finite_point_leaf.py",
            ]
        },
        verifier_sha256=digest(Path(__file__)),
        dual_sha256=digest(cert),
        cloud_sha256=digest(dst / "q-squeezed.json"),
    )
    output(out, result)
    print(guaranteed, maxloss, max(drops), flush=True)


if __name__ == "__main__":
    p = ArgumentParser(description=__doc__)
    p.add_argument("computation", type=Path)
    p.add_argument("certificate", type=Path)
    p.add_argument("out", type=Path)
    a = p.parse_args()
    run(a.computation.resolve(), a.certificate, a.out)
