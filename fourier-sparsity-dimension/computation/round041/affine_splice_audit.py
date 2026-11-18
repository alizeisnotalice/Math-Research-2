#!/usr/bin/env python3
"""Exhaustive small binary and quaternary Boolean splices with noncharacter leaves."""

import argparse
from collections import Counter
import hashlib
from itertools import product
import json
from pathlib import Path
import platform
import sys
import time


def fwht(values):
    out = list(values)
    width = 1
    while width < len(out):
        for start in range(0, len(out), 2 * width):
            for j in range(start, start + width):
                lo, hi = out[j], out[j + width]
                out[j], out[j + width] = lo + hi, lo - hi
        width *= 2
    return out


def gf2_rank(values):
    pivots = {}
    for value in values:
        while value:
            i = value.bit_length() - 1
            if i not in pivots:
                pivots[i] = value
                break
            value ^= pivots[i]
    return len(pivots)


def spectra(dim):
    points = 1 << dim
    return [
        [(x, w) for x, w in enumerate(fwht(1 - 2 * ((f >> y) & 1) for y in range(points))) if w]
        for f in range(1 << points)
    ]


def embed(frequency, shared, private, branch):
    shared_mask = (1 << shared) - 1
    return (frequency & shared_mask) | ((frequency >> shared) << (shared + private * branch))


def direct_witness(family, recorded, selector_bits, shared, private):
    """Transform the full assembled truth table, independently from leaf embeddings."""
    leaves = 1 << selector_bits
    data_dim = shared + private * leaves
    xs = range(1 << (data_dim + selector_bits))
    mask = (1 << private) - 1
    full = []
    for x in xs:
        select, data = x % leaves, x // leaves
        local = (data % (1 << shared)) | (((data >> (shared + private * select)) & mask) << shared)
        full.append(1 - 2 * ((family[select] >> local) & 1))
    walsh = fwht(full)
    denominator_scale = 1 << (data_dim - shared - private)
    assert {xi: w for xi, w in enumerate(walsh) if w} == {
        xi: w * denominator_scale for xi, w in recorded.items()
    }
    return sum(int(v < 0) << x for x, v in enumerate(full))


def audit_parallel_splice(selector_bits, shared, private):
    leaves = 1 << selector_bits
    dim = shared + private * leaves + selector_bits
    lifted = 1 << dim
    data_lifted = 1 << (dim - selector_bits)
    library = spectra(shared + private)
    leaf_spectra = [
        [[(embed(xi, shared, private, j), w) for xi, w in spec] for spec in library]
        for j in range(leaves)
    ]
    stage = Counter()
    nonchar_colliding_slacks = Counter()
    witnesses = {}

    def record(name, family, assembled, rows, ranks, terms):
        if name in witnesses:
            return
        truth_mask = direct_witness(family, assembled, selector_bits, shared, private)
        witnesses[name] = {
            "leaf_truth_masks": [hex(f) for f in family],
            "assembled_truth_mask": hex(truth_mask),
            "leaf_sparsities": [len(spec) for spec in terms],
            "colliding_frequency_rows": {
                beta: row for beta, row in rows.items() if sum(bool(x) for x in row) > 1
            },
            "assembled_spectrum": dict(sorted(assembled.items())),
            "sparsity": len(assembled),
            "lifted_rank": ranks[0],
            "projected_lifted_rank": ranks[1],
        }

    for family in product(range(len(library)), repeat=leaves):
        rows = {}
        terms = [leaf_spectra[j][f] for j, f in enumerate(family)]
        for j, spec in enumerate(terms):
            for beta, value in spec:
                rows.setdefault(beta, [0] * leaves)[j] = value
        assembled = {}
        collisions = 0
        incidence = 0
        unequal_binary = False
        tight_binary = True
        unsplit_all_equal_row = False
        row_loss = 0
        for beta, row in rows.items():
            incidence_at_beta = sum(bool(x) for x in row)
            incidence += incidence_at_beta
            collisions += incidence_at_beta > 1
            selector_spectrum = fwht(row)
            entries = [(z, w) for z, w in enumerate(selector_spectrum) if w]
            assert incidence_at_beta * len(entries) >= leaves
            row_loss += leaves - len(entries)
            if selector_bits == 1 and incidence_at_beta == 2:
                unequal_binary |= len(entries) == 2
                tight_binary &= len(entries) == 1
            if incidence_at_beta == leaves and len({abs(x) for x in row}) == 1:
                unsplit_all_equal_row |= len(entries) == leaves
            for z, w in entries:
                assembled[(beta << selector_bits) | z] = w
        actual = len(assembled)
        data_rank = gf2_rank(beta | data_lifted for beta in rows)
        combined_rank = gf2_rank(xi | lifted for xi in assembled)
        selector_rank = combined_rank - data_rank
        assert 0 <= selector_rank <= selector_bits
        assert actual == leaves * len(rows) - row_loss
        assert actual * incidence >= leaves * len(rows) ** 2
        if selector_bits == 1:
            al, ar = [gf2_rank(beta | data_lifted for beta, _ in spec) for spec in terms]
            h = al + ar - data_rank
            assert collisions <= (1 << (h - 1) if h else 0)
            naive = 2 * incidence
            assert 2 * collisions <= naive - actual <= 3 * collisions
            if any(len(spec) == 1 for spec in terms) and any(len(spec) == 4 for spec in terms):
                assert 4 * actual >= (combined_rank + 1) ** 2

        ranks = combined_rank, data_rank
        stage["families"] += 1
        any_noncharacter = any(len(spec) > 1 for spec in terms)
        all_noncharacter = all(len(spec) > 1 for spec in terms)
        stage["with_noncharacter"] += any_noncharacter
        stage["all_noncharacter"] += all_noncharacter
        stage["with_collisions"] += collisions > 0
        if not (collisions and any_noncharacter):
            continue

        score = 4 * actual - (combined_rank + 1) ** 2
        nonchar_colliding_slacks[score] += 1
        stage["noncharacter_and_colliding"] += 1
        stage["all_noncharacter_and_colliding"] += all_noncharacter
        if score < 0:
            record("violates_quadratic_budget", family, assembled, rows, ranks, terms)
        if score == 0 and all_noncharacter:
            record("tight_with_all_noncharacter_and_collision", family, assembled, rows, ranks, terms)
        if selector_bits == 1:
            if all_noncharacter and collisions ** 2 > min(len(spec) for spec in terms):
                record("collision_number_exceeds_smaller_square_root", family, assembled, rows, ranks, terms)
            if all_noncharacter and unequal_binary:
                record("binary_collision_without_coefficient_cancellation", family, assembled, rows, ranks, terms)
            if all_noncharacter and tight_binary and len(rows) == collisions:
                if selector_rank == 0 and family[0] != family[1]:
                    record("single_fiber_per_frequency_affine_graph", family, assembled, rows, ranks, terms)
                if selector_rank == 1:
                    record("single_fiber_per_frequency_nonaffine_graph", family, assembled, rows, ranks, terms)
        if selector_bits == 2 and all_noncharacter and unsplit_all_equal_row:
            record("quaternary_equal_magnitude_collision_without_cancellation", family, assembled, rows, ranks, terms)

    return {
        "selector_bits": selector_bits,
        "shared_data_bits": shared,
        "private_bits_per_leaf": private,
        "ambient_bits": dim,
        "counts": dict(stage),
        "quadratic_slack_histogram_for_colliding_noncharacter": dict(sorted(nonchar_colliding_slacks.items())),
        "witnesses": witnesses,
    }


def audit_nonparallel():
    points = range(16)
    blocks = [
        [x for x in points if x & 1 == b and (x >> (1 if b == 0 else 2)) & 1 == a]
        for b, a in product((0, 1), repeat=2)
    ]
    assert sorted(x for block in blocks for x in block) == list(points)
    signs = [sum(((f >> i) & 1) << x for i, x in enumerate(block)) for block in blocks for f in range(16)]
    branches = [[fwht((-1) ** ((f >> i) & 1) if x in block else 0 for x in points
                      for i in [block.index(x) if x in block else 0]) for f in range(16)] for block in blocks]
    result = Counter()
    witness = None
    for family in product(range(16), repeat=4):
        noncharacter = [f.bit_count() % 2 == 1 for f in family]
        if not any(noncharacter):
            continue
        summands = [branches[j][f] for j, f in enumerate(family)]
        w = [sum(col) for col in zip(*summands)]
        support = [x for x in points if w[x]]
        rank = gf2_rank(x | 16 for x in support)
        collisions = sum(sum(bool(value) for value in col) > 1 for col in zip(*summands))
        cancelled_collisions = sum(w[x] == 0 and sum(bool(v[x]) for v in summands) > 1 for x in points)
        score = 4 * len(support) - (rank + 1) ** 2
        result[(sum(noncharacter), score, collisions, cancelled_collisions)] += 1
        if witness is None and all(noncharacter) and score == 0 and collisions:
            truth_mask = sum(signs[16 * j + f] for j, f in enumerate(family))
            assert w == fwht(1 - 2 * ((truth_mask >> x) & 1) for x in points)
            witness = {
                "leaf_truth_masks": [hex(f) for f in family],
                "assembled_truth_mask": hex(truth_mask),
                "assembled_spectrum": w,
                "collision_count": collisions,
                "cancelled_collision_count": cancelled_collisions,
            }
    return {
        "blocks": blocks,
        "details_by_noncharacter_count_slack_collisions_cancelled": [
            [*key, value] for key, value in sorted(result.items())
        ],
        "tight_witness": witness,
        "noncharacter_splices": sum(result.values()),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    start = time.perf_counter()
    shapes = [(1, i, 3 - i) for i in range(4)] + [(2, i, 2 - i) for i in range(3)]
    result = {
        "python": platform.python_version(),
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "arithmetic": "Exact integers; deterministic exhaustive enumeration.",
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "parallel_affine_splices": [audit_parallel_splice(*shape) for shape in shapes],
        "nonparallel_affine_partition": audit_nonparallel(),
        "elapsed_seconds": time.perf_counter() - start,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print("Checked all seven parallel and one nonparallel partition families.")


if __name__ == "__main__":
    main()
