#!/usr/bin/env python3
"""Rank and density census for complete affine line matchings in F_2^3."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from functools import cache
import hashlib
from itertools import combinations
import json
from pathlib import Path
import platform
import time


N = 3
ADDR = 1 << N


def members(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def word_mask(points):
    return sum(1 << x for x in points)


def parity(x):
    return x.bit_count() & 1


def rank(words):
    basis = {}
    for v in words:
        while v:
            p = v.bit_length() - 1
            if p not in basis:
                basis[p] = v
                break
            v ^= basis[p]
    return len(basis)


def span(generators):
    result = {0}
    for g in generators:
        result |= {x ^ g for x in result}
    return result


def gaussian(n, k):
    if not 0 <= k <= n:
        return 0
    q = [1] + [0] * k
    for i in range(1, n + 1):
        for j in range(min(i, k), 0, -1):
            q[j] += (1 << (i - j)) * q[j - 1]
    return q[k]


@cache
def rref_subspaces(n):
    space_masks = []
    counts = Counter()
    for r in range(n + 1):
        for pivot in combinations(range(n), r):
            free = [(i, j) for i, p in enumerate(pivot) for j in range(p + 1, n) if j not in pivot]
            for k in range(1 << len(free)):
                rows = [1 << p for p in pivot]
                for bit, (i, j) in enumerate(free):
                    rows[i] |= ((k >> bit) & 1) << j
                W = span(rows)
                assert len(W) == 1 << r
                counts[r] += 1
                space_masks.append((r, word_mask(W)))
    assert len({w for _, w in space_masks}) == len(space_masks)
    assert all(counts[j] == gaussian(n, j) for j in range(n + 1))
    return tuple(space_masks)


def matchings(points):
    if not points:
        yield ()
        return
    first, *tail = points
    for j, mate in enumerate(tail):
        for rest in matchings(tail[:j] + tail[j + 1:]):
            yield ((first, mate),) + rest


def characters_partition(blocks):
    """Compute the complete actual Walsh table of the selector with private characters."""
    lookup = {a: i for i, block in enumerate(blocks) for a in block}
    vals = [1 - 2 * ((x >> (N + lookup[x % ADDR])) & 1) for x in range(ADDR << len(blocks))]
    w = vals[:]
    step = 1
    while step < len(w):
        for left in range(0, len(w), 2 * step):
            for i in range(left, left + step):
                a, b = w[i], w[i + step]
                w[i], w[i + step] = a + b, a - b
        step *= 2
    return w


def geometry(blocks):
    assert sorted(x for block in blocks for x in block) == list(range(ADDR))
    normal_spaces = []
    for block in blocks:
        origin = min(block)
        direction = {x ^ origin for x in block}
        assert direction == span(direction)
        normal = {alpha for alpha in range(ADDR) if all(not parity(alpha & v) for v in direction)}
        assert len(normal) * len(block) == ADDR
        normal_spaces.append(normal)
    m = len(blocks)
    normal_generators = set().union(*normal_spaces)
    W = span(normal_generators)
    R = rank(normal_generators)
    assert len(W) == 1 << R
    wmask = word_mask(W)
    cells = [word_mask((1 << (N + i)) ^ a for a in ws) for i, ws in enumerate(normal_spaces)]
    true_walsh = characters_partition(blocks)
    support = word_mask(i for i, w in enumerate(true_walsh) if w)
    assert support == sum(cells)
    assert rank(members(support)) == R + m
    common_scale = 1 << m
    for i, cell in enumerate(cells):
        assert all(abs(true_walsh[x]) == len(blocks[i]) * common_scale for x in members(cell))
    return normal_spaces, R, wmask, cells, support


def explicit_profiles(blocks, spaces, R, wmask, *, strict):
    m = len(blocks)
    delta = m - R
    assert delta >= 1
    totals = Counter()
    witness = {}
    table = {}
    subspaces = [(d, q) for d, q in rref_subspaces(N) if q & wmask == q]
    assert len(subspaces) == sum(gaussian(R, k) for k in range(R + 1))
    normal_points = [word_mask(w) for w in spaces]
    for mask in range((1 << m) - 1):
        group = tuple(members(mask))
        actual = rank(a for i in group for a in spaces[i])
        deficiency = len(group) - actual
        totals["proper_subfamily_deficiency_checks"] += 1
        totals["proper_subfamily_deficiency_equalities"] += deficiency == delta - 1
        if deficiency >= delta:
            totals["proper_subfamily_deficiency_violations"] += 1
            witness.setdefault("proper_deficiency_violation", {
                "J": group, "deficiency": deficiency, "delta": delta,
            })
        if strict:
            assert deficiency <= delta - 1
    augmented_bit = 1 << (N + m)
    for qdim, q in subspaces:
        fibres = [tuple(members(w & q)) for w in normal_points]
        for subset in range(1 << m):
            indices = tuple(members(subset))
            points = tuple((1 << (N + i)) ^ a for i in indices for a in fibres[i])
            r = rank(a ^ augmented_bit for a in points)
            assert r == rank(points)
            assert r == len(indices) + rank(a for i in indices for a in fibres[i])
            assert r <= qdim + len(indices)
            n = len(points)
            scaled_weight = sum(len(blocks[i]) * len(fibres[i]) for i in indices)
            table[q, subset] = n, r, scaled_weight
            totals["Q_and_J_profiles"] += 1
            totals["Q_and_J_with_rank_below_upper_budget"] += r < qdim + len(indices)
            if r == m + R:
                totals["full_rank_profiles"] += 1
                assert n == len(indices) * (ADDR // 2) if strict else scaled_weight == ADDR * m
                continue
            rhs = ADDR * (r + delta - 1)
            residual = rhs - 2 * scaled_weight
            if strict:
                assert 2 * n <= m * (r + delta - 1)
                assert residual == (ADDR // m) * (m * (r + delta - 1) - 2 * n)
            if residual < 0:
                totals["weighted_density_violations"] += 1
                witness.setdefault("weighted_density_violation", {
                    "Q": list(members(q)), "J": indices, "rank": r,
                    "weighted_numerator_out_of_8": scaled_weight, "delta": delta, "density_residual": residual,
                })
            if residual == 0:
                totals["proper_weighted_density_equalities"] += 1
                witness.setdefault("proper_weighted_density_equality", {
                    "Q": list(members(q)), "J": indices, "rank": r,
                    "weighted_numerator_out_of_8": scaled_weight, "delta": delta,
                })
    return table, totals, witness


def verify_full_ambient(space_data, wmask, cells, support, table, blocks, R):
    m = len(cells)
    augmented_bit = 1 << (N + m)
    counts = Counter()
    for ambient_r, V in space_data:
        q = V & wmask
        part = V & support
        fragments = [part & cell for cell in cells]
        J = sum((1 << i) for i, restricted in enumerate(fragments) if restricted)
        mass, r, weighted_mass = table[q, J]
        assert part.bit_count() == mass
        assert rank(x ^ augmented_bit for x in members(part)) == r
        assert weighted_mass == sum(fragment.bit_count() * len(blocks[i]) for i, fragment in enumerate(fragments))
        assert r <= ambient_r
        normals = []
        reps = []
        for i, fragment in enumerate(fragments):
            if fragment:
                origin = (fragment & -fragment).bit_length() - 1
                reps.append(origin)
                projected_normal_points = {x ^ origin for x in members(fragment)}
                assert projected_normal_points == set(members(q)).intersection(
                    x ^ (1 << (N + i)) for x in members(cells[i])
                )
                normals.extend(projected_normal_points)
        assert rank(reps + normals) == r
        counts["actual_ambient_subspaces"] += 1
        counts["actual_proper_rank_subspaces"] += r < m + R
        counts["actual_full_rank_subspaces"] += r == m + R
    return counts


def matching_job(blocks):
    ws, R, wmask, cells, support = geometry(blocks)
    table, totals, witness = explicit_profiles(blocks, ws, R, wmask, strict=True)
    totals.update(verify_full_ambient(rref_subspaces(N + 4), wmask, cells, support, table, blocks, R))
    totals["matching_partitions"] = 1
    return R, tuple(sorted(Counter(a ^ b for a, b in blocks).values())), dict(totals), witness


def unequal_cases():
    shapes = [
        ((0, 1, 2, 3), (4, 5), (6, 7)),
        ((0, 1, 2, 3), (4, 6), (5,), (7,)),
        ((0, 1, 2, 3), (4,), (5,), (6,), (7,)),
        ((0, 1), (2, 4), (3, 7), (5,), (6,)),
        ((0, 1), (2, 4), (3,), (5,), (6,), (7,)),
        ((0, 1), (2,), (3,), (4,), (5,), (6,), (7,)),
    ]
    reports = []
    for blocks in shapes:
        ws, R, wmask, cells, support = geometry(blocks)
        table, totals, witness = explicit_profiles(blocks, ws, R, wmask, strict=False)
        reports.append({
            "blocks": blocks, "normal_rank": R, "delta": len(blocks) - R,
            "counts": dict(totals), "witness": witness,
        })
    return reports


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--workers", type=int, default=12)
    args = parser.parse_args()
    started = time.perf_counter()
    all_matchings = list(matchings(list(range(ADDR))))
    assert len(all_matchings) == 105
    totals = Counter()
    group = {}
    witnesses = {}
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for R, multiplicity, vals, witness in pool.map(matching_job, all_matchings, chunksize=2):
            key = (R, multiplicity)
            group.setdefault(key, Counter()).update(vals)
            totals.update(vals)
            if key not in witnesses:
                witnesses[key] = witness
    unequal = unequal_cases()
    report = {
        "python": platform.python_version(), "platform": platform.platform(),
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "workers": args.workers, "totals": dict(totals),
        "matching_groups_by_R_and_direction_multiplicity": [
            {"R": R, "delta": 4 - R, "direction_multiplicity": multi, "counts": dict(group[R, multi]),
             "witness": witnesses[R, multi]}
            for R, multi in sorted(group)
        ],
        "weighted_unequal_affine_partitions": unequal,
        "elapsed_seconds": time.perf_counter() - started,
    }
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report["totals"], sort_keys=True))


if __name__ == "__main__":
    main()
