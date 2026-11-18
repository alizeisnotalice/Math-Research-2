"""Certify global area growth under arbitrary parallel translations of four bodies."""

from argparse import ArgumentParser
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import platform

import flint
from flint import arb, arb_mat, ctx
from certify_five_body_contact_stationary import (
    BOX, add, det, interval, left, model, normal, q, scale, sub,
)
from certify_five_body_contact_transverse import EVENT_ORDER, SIMPLEX, face_orders

HERE = Path(__file__).resolve().parent
WEIGHTS = HERE / "five-body-parallel-weights.json"
SMALL = q(".0009661")
DEN = 10**12


def velocities(owner, point, pose):
    speeds = [(0, 0)] * 9
    if owner:
        j = owner - 3
        speeds[j:j + 2] = [(1, 0), (0, 1)]
        if owner == 9:
            speeds[8] = left(sub(point, pose[9:11]))
    return speeds


def connector(u, X, w, Y, pose):
    speeds = zip(velocities(u, X, pose), velocities(w, Y, pose))
    return [det(add(a, b), sub(Y, X)) / 2 for a, b in speeds]


def chain_gradients():
    bounds, pose, geometry, _ = model(*(interval(*row) for row in BOX))
    normals = {tag: normal(end) for tag, end, *_ in bounds if tag in EVENT_ORDER}
    orders = face_orders(geometry, normals)
    base = [arb(0)] * 9
    increments = {}
    start = arb(0)
    for piece, following in zip(bounds, bounds[1:] + bounds[:1]):
        tag, end, u, _, rad, center = piece
        _, _, w, _, next_rad, next_center = following
        n = normal(end)
        diff = sub(n, normal(start))
        for j, v in enumerate(velocities(u, center, pose)):
            base[j] += rad * det(v, diff)
        X = add(center, scale(rad, n))
        Y = add(next_center, scale(next_rad, n))
        original = connector(u, X, w, Y, pose)
        base = [a + b for a, b in zip(base, original)]
        if tag in orders:
            seq = orders[tag]
            assert seq[0][0] == u and seq[-1][0] == w
            inc = []
            for flags in product((0, 1), repeat=len(seq) - 2):
                path = [seq[0], *(p for p, f in zip(seq[1:-1], flags) if f), seq[-1]]
                pieces = [connector(*a, *b, pose) for a, b in zip(path, path[1:])]
                inc.append([sum(row[j] for row in pieces) - original[j] for j in range(9)])
            increments[tag] = inc
        start = end
    return base, [increments[k] for k in EVENT_ORDER]


def small_cube(base, increments):
    gradients = [
        [base[j] + sum(event[k][j] for event, k in zip(increments, picks)) for j in range(9)]
        for picks in SIMPLEX
    ]
    matrix = arb_mat([[1] * 10, *[[g[j] for g in gradients] for j in range(9)]])
    assert abs(matrix.det()) > 0
    inverse = matrix.inv()
    for i in range(10):
        assert inverse[i, 0] > q(".05")
        assert inverse[i, 0] > SMALL * sum(abs(inverse[i, j]) for j in range(1, 9))
    return str(matrix.det()), [str(inverse[i, 0]) for i in range(10)]


def exact_weights(values, cuts):
    groups = []
    for low, high in zip(cuts, cuts[1:]):
        row = values[low:high]
        weights = [max(0, round(DEN * v)) for v in row]
        largest = max(range(len(row)), key=lambda i: row[i])
        weights[largest] += DEN - sum(weights)
        assert min(weights) >= 0
        groups.append([[j, weight] for j, weight in enumerate(weights) if weight])
    return groups


def projected_gradient(base, increments, groups):
    y = list(base[:8])
    assert len(groups) == len(increments)
    for choices, group in zip(increments, groups):
        assert sum(w for k, w in group) == DEN
        for k, w in group:
            assert 0 <= k < len(choices) and 0 < w <= DEN
            for j in range(8):
                y[j] += choices[k][j] * w / DEN
    return y


def locate(base, increments):
    import numpy as np
    from scipy.optimize import linprog

    lengths = [len(e) for e in increments]
    cuts = np.cumsum([0, *lengths])
    m = cuts[-1]
    A = np.zeros((14, m + 1))
    A[6:, :m] = [[float(g[j].mid()) for e in increments for g in e] for j in range(8)]
    for i, (l, h) in enumerate(zip(cuts, cuts[1:])):
        A[i, l:h] = 1
    rhs = np.r_[np.ones(6), [-float(v.mid()) for v in base[:8]]]
    obj = np.r_[np.zeros(m), -1]
    witnesses = []
    for signs in product((-1, 1), repeat=8):
        A[6:, -1] = -np.asarray(signs)
        res = linprog(obj, A_eq=A, b_eq=rhs, bounds=(0, None), method="highs")
        assert res.success
        witnesses.append([round(DEN * res.x[-1]), exact_weights(res.x, cuts)])
    # Leave sufficient slack for the interval correction via the seed cube.
    target = min(w[0] for w in witnesses) // 10**6 - 1
    data = dict(denominator=DEN, target_millionths=int(target), rows=witnesses)
    # For common-direction motions it suffices to majorize four signed
    # projections. A single witness covers both endpoints and hence every
    # direction between them by linear interpolation.
    rank_cases = {}
    g = A[6:, :m].reshape((4, 2, m))
    shifts = -rhs[6:].reshape((4, 2))
    for vertical in (-1, 1):
        for signs in product((-1, 1), repeat=4):
            ss = np.asarray(signs)
            stack = [(0, 1)]
            rows = []
            while stack:
                k, n = stack.pop()
                ends = np.asarray([(v / n, vertical * (1 - v / n)) for v in (k, k + 1)])
                projected = np.einsum("ej,ijk->eik", ends, g)
                fixed = np.einsum("ej,ij->ei", ends, shifts)
                ub = np.c_[(-projected * ss[None, :, None]).reshape((8, m)), np.ones(8)]
                res = linprog(obj, A_ub=ub, b_ub=(fixed * ss).ravel(), A_eq=A[:6],
                              b_eq=np.ones(6), bounds=(0, None), method="highs")
                assert res.success
                if res.x[-1] < .00240001:
                    assert n < 256
                    stack.extend([(2*k+1, 2*n), (2*k, 2*n)])
                    continue
                rows.append([k, n, exact_weights(res.x, cuts)])
            rank_cases[str((vertical, *signs))] = rows
    data["parallel_direction_cases"] = rank_cases
    WEIGHTS.write_text(json.dumps(data, separators=(",", ":")) + "\n")
    print("certification target", target, "/ 1000000")


def certify(base, inc):
    det_ball, barycenters = small_cube(base, inc)
    data = json.loads(WEIGHTS.read_bytes())
    R = Fraction(data["target_millionths"], 10**6)
    target = q(str(R))
    assert target > SMALL and data["denominator"] == DEN and len(data["rows"]) == 256
    max_error = arb(0)
    for signs, (r_num, groups) in zip(product((-1, 1), repeat=8), data["rows"]):
        radius = q(str(Fraction(r_num, DEN)))
        assert radius > target and len(groups) == len(inc)
        y = projected_gradient(base, inc, groups)
        # The exact convex combination (R/r)y + (1-R/r)z covers this vertex
        # of the large cube for some z in the certified seed cube.
        for j, s in enumerate(signs):
            err = abs(radius * s - y[j])
            assert target * err < (radius - target) * SMALL
            max_error = max(max_error, err.upper())
    parallel = q(".0024")
    rank_cases = data["parallel_direction_cases"]
    parallel_count = 0
    min_margin = None
    for vertical in (-1, 1):
        for signs in product((-1, 1), repeat=4):
            cases = rank_cases[str((vertical, *signs))]
            next_start = Fraction(0)
            for k, n, groups in cases:
                assert n > 0 and Fraction(k, n) == next_start and len(groups) == len(inc)
                next_start = Fraction(k + 1, n)
                y = projected_gradient(base, inc, groups)
                for x in (k, k + 1):
                    for i, sign in enumerate(signs):
                        proj = sign * (x*y[2*i] + vertical*(n-x)*y[2*i+1]) / n
                        assert proj > parallel
                        margin = (proj-parallel).lower()
                        min_margin = margin if min_margin is None else min(min_margin, margin)
                parallel_count += 1
            assert next_start == 1
    threshold = q(".837295006222")
    width = q("1.245477")
    width_area = (arb.pi() / 4 + ((width * width - 1).sqrt() - (1 / width).acos()) / 2)
    assert width_area > threshold
    rho = q("1e-6") * (arb.pi() / 10).tan() / 2
    delta = q("2e-7")
    example_gain = parallel / 1000 - (width + 2 * rho) * delta - 4 * threshold * rho - arb.pi() * rho * rho
    assert example_gain > q("1.6e-6")
    paths = [Path(__file__), WEIGHTS, HERE / "certify_five_body_contact_stationary.py",
             HERE / "certify_five_body_contact_transverse.py",
             HERE / "certify_four_body_contact_stationary.py"]
    print(json.dumps(dict(
        python=platform.python_version(), flint=flint.__version__, precision_bits=ctx.prec,
        source_sha256={p.name: sha256(p.read_bytes()).hexdigest() for p in paths},
        normal_box=BOX, directions=["t3x", "t3y", "t5x", "t5y", "t7x", "t7y", "t9x", "t9y"],
        auxiliary_direction="theta9", seed_simplex=SIMPLEX, simplex_determinant=det_ball,
        seed_barycentric=barycenters, seed_cube_radius="9661/10000000",
        cube_radius=str(R), vertices=256, maximum_vertex_rounding_error=str(max_error),
        global_parallel_coefficient="3/1250", angular_cells=parallel_count,
        minimum_projection_slack=str(min_margin), width_at_threshold="1.245477", minimum_area_at_width=str(width_area),
        example_contradiction_gain=str(example_gain),
        scope="All common-direction translations at the five-body contact stationary point",
    ), indent=2))


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--locate", action="store_true")
    args = parser.parse_args()
    base, inc = chain_gradients()
    if args.locate:
        locate(base, inc)
    else:
        certify(base, inc)
