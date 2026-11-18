#!/usr/bin/env python3
"""Integer Walsh, weighted frequency incidence, and projected rank checks."""

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from functools import lru_cache
import hashlib
from itertools import permutations, product
import json
from math import lcm
from pathlib import Path
import os
import platform
import time


ADDR = 8
DATA_BITS = 4
DATA = 1 << DATA_BITS
AMBIENT = ADDR * DATA
SCALE = lcm(*range(1, ADDR + 1))
SQ_SCALE = SCALE * SCALE


def parity(n):
    return n.bit_count() & 1


def rank(words):
    root = {}
    for x in words:
        while x:
            p = x.bit_length() - 1
            if p not in root:
                root[p] = x
                break
            x ^= root[p]
    return len(root)


def fwht(values):
    out = list(values)
    width = 1
    while width < len(out):
        for i in range(0, len(out), 2 * width):
            for j in range(i, i + width):
                a, b = out[j], out[j + width]
                out[j], out[j + width] = a + b, a - b
        width *= 2
    return out


def affine_coordinates(points):
    origin = min(points)
    ordered = [origin]
    basis = []
    for point in sorted(points):
        if point not in ordered:
            v = point ^ origin
            basis.append(v)
            ordered.extend(p ^ v for p in ordered[:])
    assert set(ordered) == set(points)
    return ordered, basis


def affine_set(points):
    origin = min(points)
    return len(points) == 1 << rank(x ^ origin for x in points)


def phase_choices(d):
    def e(i):
        return 1 << (d + i)

    local = 1 if d else e(3)
    return [
        (0, 0, 0, 0), (0, 0, e(0), 0), (0, 0, e(1), 1),
        (0, 0, e(0) ^ e(1), 0), (0, 0, e(2), 0),
        (local, e(0), e(1), 0), (local, e(1), e(2), 0), (local, e(2), e(0), 1),
        (e(0), e(1), e(2), 0), (e(0), e(2), 0, 0), (e(1), e(2), e(0), 1),
        (local ^ e(0), e(1), e(3), 0),
    ]


class Partition:
    def __init__(self, blocks):
        assert sorted(x for block in blocks for x in block) == list(range(ADDR))
        self.blocks = blocks
        self.coordinates = [affine_coordinates(set(block)) for block in blocks]
        direction_sets = [{x ^ min(block) for x in block} for block in blocks]
        common = set.intersection(*direction_sets)
        self.shared_order, self.address_basis = affine_coordinates(common)
        self.projection = [
            self.project(x, self.address_basis, ADDR.bit_length() - 1) for x in range(AMBIENT)
        ]
        self.width = len(self.address_basis) + DATA_BITS
        self.lift = 1 << self.width
        self.affine_cache = {}

    def project(self, f, forms, local_address_bits):
        lower = sum(parity(f & v) << j for j, v in enumerate(forms))
        return lower | ((f >> local_address_bits) << len(forms))

    def make_leaf(self, i, phase):
        coordinates, local_basis = self.coordinates[i]
        d = len(local_basis)
        q, r, linear, constant = phase
        if q or r:
            assert q and r and q != r
        local_sign = [
            1 - 2 * ((parity(q & x) & parity(r & x)) ^ parity(linear & x) ^ constant)
            for x in range(DATA << d)
        ]
        full_mask = sum(
            int(value < 0) << (coordinates[x % len(coordinates)] | ((x >> d) << 3))
            for x, value in enumerate(local_sign)
        )
        zero_local = coordinates[0]
        in_local = [coordinates.index(x ^ zero_local) for x in self.address_basis]
        w = fwht(local_sign)
        support = tuple(x for x, value in enumerate(w) if value)
        assert len(support) == (4 if q else 1)
        if q:
            assert rank(x ^ support[0] for x in support) == 2
        projection = Counter(self.project(x, in_local, d) for x in support)
        return full_mask, projection, support, local_sign

    def active_union_affine(self, active):
        key = tuple(active)
        if key not in self.affine_cache:
            self.affine_cache[key] = affine_set(set().union(*(self.blocks[i] for i in key)))
        return self.affine_cache[key]


@lru_cache(maxsize=None)
def weighted_basis(lifted):
    best = []
    for vector, weight in sorted(lifted, key=lambda r: (-r[1], r[0])):
        if rank(best + [vector]) > len(best):
            best.append(vector)
    weights = dict(lifted)
    greedy = sum(weights[v] for v in best)
    optimum = 0
    for flags in range(1 << len(lifted)):
        members = [r for i, r in enumerate(lifted) if (flags >> i) & 1]
        if rank(v for v, _ in members) == len(members):
            optimum = max(optimum, sum(w for _, w in members))
    assert greedy == optimum
    return greedy, tuple(best)


def audit(geo, leaves, *, detailed=False):
    total_truth = 0
    beta_incidence = {}
    for i, (mask, proj, _, _) in enumerate(leaves):
        total_truth |= mask
        for beta, count in proj.items():
            beta_incidence.setdefault(beta, []).append((i, count))
    spectrum = fwht(1 - 2 * ((total_truth >> x) & 1) for x in range(AMBIENT))
    support = [x for x, val in enumerate(spectrum) if val]
    s = len(support)
    frequencies = Counter(geo.projection[x] for x in support)
    assert set(frequencies) == set(beta_incidence)
    measures = [len(a) for a in geo.blocks]
    budget = [0] * len(leaves)
    rank_item = [dict() for part in leaves]
    shared_cases = private_cases = multi_private = multi_shared = loose_private_formula = 0
    failures, maximum_failures, shared_affine_count = 0, 0, 0
    record = {}
    for beta, incidences in beta_incidence.items():
        mass = sum(measures[i] for i, _ in incidences)
        sb = frequencies[beta]
        assert sb * mass >= ADDR
        if len(incidences) == 1:
            i, r = incidences[0]
            assert sb * measures[i] == r * ADDR
            budget[i] += r * ADDR * SQ_SCALE // measures[i]
            private_cases += 1
            multi_private += r > 1
            loose_private_formula += r > 1 and sb * measures[i] > ADDR
            rank_item[i][beta | geo.lift] = 0
        else:
            real_r_sum = sum(measures[i] * r for i, r in incidences)
            shared_cases += 1
            multi_shared += any(r > 1 for _, r in incidences)
            affine = geo.active_union_affine(i for i, _ in incidences)
            local_max = max(r for _, r in incidences)
            if affine:
                assert sb * mass >= ADDR * local_max
                shared_affine_count += 1
            if sb * mass * mass < ADDR * real_r_sum:
                assert not affine
                failures += 1
            if sb * mass < ADDR * local_max:
                assert not affine
                maximum_failures += 1
            for i, _ in incidences:
                budget[i] += ADDR * measures[i] * SQ_SCALE // (mass * mass)
                rank_item[i][beta | geo.lift] = SCALE * (mass - measures[i]) // mass
        if detailed:
            partial = fwht(val if geo.projection[x] == beta else 0 for x, val in enumerate(spectrum))
            occupied = [x for x, val in enumerate(partial) if val]
            assert occupied and len(occupied) * sb >= AMBIENT
            active_addresses = set().union(*(geo.blocks[i] for i, _ in incidences))
            assert all(x % ADDR in active_addresses for x in occupied)
            record[beta] = {
                "band_support": sb, "P": str(Fraction(mass, ADDR)),
                "actual_support_probability": str(Fraction(len(occupied), AMBIENT)),
                "participants_local_preimage_counts": incidences,
            }
    real_budget = sum(budget)
    assert real_budget <= SQ_SCALE * s
    deduct = [weighted_basis(tuple(sorted(w.items()))) for w in rank_item]
    local_rank = [rank(w.keys()) for w in rank_item]
    union_rank = rank(beta | geo.lift for beta in frequencies)
    rank_bound = sum(SCALE * d - b for d, (b, _) in zip(local_rank, deduct))
    assert SCALE * union_rank <= rank_bound
    if not detailed:
        return {
            "cases": 1, "walsh_coefficients": AMBIENT,
            "bands_checked": len(frequencies), "private_bands": private_cases,
            "private_bands_with_multiple_preimages": multi_private,
            "private_bands_where_only_1_over_P_is_strict": loose_private_formula,
            "shared_bands": shared_cases, "shared_bands_with_multiple_preimages": multi_shared,
            "shared_bands_affine_union": shared_affine_count,
            "invalid_unconditional_weighted_multiplicity_band_failures": failures,
            "invalid_unconditional_maximal_multiplicity_band_failures": maximum_failures,
            "exact_total_sparsity_budget": real_budget == SQ_SCALE * s,
            "exact_weighted_projection_rank_budget": SCALE * union_rank == rank_bound,
        }
    expected = [Fraction(0)] * len(leaves)
    probability = Fraction(0)
    for order in permutations(range(len(leaves))):
        prev = set()
        remain = ADDR
        p = Fraction(1)
        accrued = []
        for i in order:
            p *= Fraction(measures[i], remain)
            remain -= measures[i]
            inherited = prev.intersection(rank_item[i])
            accrued.append((i, rank(inherited)))
            prev.update(rank_item[i])
        probability += p
        for i, value in accrued:
            expected[i] += p * value
    assert probability == 1
    assert all(x >= Fraction(dd[0], SCALE) for x, dd in zip(expected, deduct))
    assert Fraction(union_rank) <= sum(Fraction(a) - x for a, x in zip(local_rank, expected))
    return {
        "blocks": geo.blocks, "intersection_address_basis": geo.address_basis,
        "truth_mask": hex(total_truth), "sparsity": s, "support": support,
        "band_details": dict(sorted(record.items())), "r_i": local_rank,
        "b_i": [str(Fraction(x, SCALE)) for x, _ in deduct],
        "chosen_weighted_independent_sets_lifted": [list(vec) for _, vec in deduct],
        "exact_expected_inherited_order_ranks": [str(x) for x in expected],
        "v_i": [str(Fraction(c * m, SQ_SCALE * ADDR)) for c, m in zip(budget, measures)],
        "sparsity_lower_bound": str(Fraction(real_budget, SQ_SCALE)),
        "projected_lifted_rank": union_rank, "projected_lifted_rank_upper_bound": str(Fraction(rank_bound, SCALE)),
        "unconditional_weighted_multiplicity_failures": failures,
        "local_spectra": [p[2] for p in leaves],
        "walsh_coefficients_at_support": [spectrum[i] for i in support],
    }


def matching(points):
    if not points:
        yield ()
        return
    x, *rest = points
    for i, y in enumerate(rest):
        for tail in matching(rest[:i] + rest[i + 1:]):
            yield ((x, y),) + tail


def task(job):
    name, blocks, indexes, first = job
    begin = time.process_time()
    geo = Partition(blocks)
    libraries = [
        [geo.make_leaf(i, phase_choices(len(basis))[j]) for j in indexes]
        for i, (_, basis) in enumerate(geo.coordinates)
    ]
    if first is not None:
        libraries[0] = [libraries[0][first]]
    totals = Counter()
    for phases in product(*libraries):
        totals.update(audit(geo, phases))
    return name, os.getpid(), time.process_time() - begin, dict(totals)


def explicit_examples():
    interlaced = ((0, 1), (2, 3), (4, 6), (5, 7))
    obj = Partition(interlaced)
    d = 1
    e = lambda j: 1 << (d + j)
    leaves = [
        (1, e(0) ^ e(1), e(0), 0),
        (1, e(1) ^ e(2), e(1), 0),
        (1, e(2) ^ e(3), e(2), 0),
        (0, 0, e(0), 0),
    ]
    tight = audit(obj, [obj.make_leaf(i, phase) for i, phase in enumerate(leaves)], detailed=True)
    assert tight["projected_lifted_rank"] == 4
    assert tight["projected_lifted_rank_upper_bound"] == "4"
    mixed_choices = [
        (1, e(0), e(1), 0),
        (0, 0, e(1), 0),
        (e(0), e(2), e(1), 0),
        (1, e(2), e(0), 0),
    ]
    mixed = audit(obj, [obj.make_leaf(i, phase) for i, phase in enumerate(mixed_choices)], detailed=True)
    assert mixed["band_details"][2]["P"] == "3/4"

    nonuniform = ((0, 1, 2, 3), (4, 6), (5,), (7,))
    geo = Partition(nonuniform)
    recipes = [(1, 4, 0, 0), (0, 0, 0, 0), (0, 0, 1, 0), (0, 0, 1, 0)]
    boundary = audit(geo, [geo.make_leaf(i, phase) for i, phase in enumerate(recipes)], detailed=True)
    assert boundary["sparsity"] == 4 and boundary["unconditional_weighted_multiplicity_failures"] == 2
    assert all(val["band_support"] == 2 and val["P"] == "3/4" for val in boundary["band_details"].values())
    output = {
        "interlaced_noninjective_tight_projected_rank": tight,
        "interlaced_characters_and_different_four_point_planes": mixed,
        "non_affine_active_union_breaks_multiplicity_weighting": boundary,
    }
    for row, recipe in zip(output.values(), (leaves, mixed_choices, recipes)):
        row["leaf_phases_q_r_linear_constant"] = recipe
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--workers", type=int, default=12)
    args = parser.parse_args()
    start = time.perf_counter()
    geometries = {
        "interlaced_two_directions": ((0, 1), (2, 3), (4, 6), (5, 7)),
        "four_distinct_line_directions": ((0, 1), (2, 4), (3, 7), (5, 6)),
        "parallel_lines_common_address": ((0, 1), (2, 3), (4, 5), (6, 7)),
        "unequal_plane_line_and_points": ((0, 1, 2, 3), (4, 6), (5,), (7,)),
        "unequal_plane_and_parallel_lines": ((0, 1, 2, 3), (4, 5), (6, 7)),
    }
    tasks = [(name, blocks, range(12), first) for name, blocks in geometries.items() for first in range(12)]
    all_matchings = list(matching(list(range(ADDR))))
    assert len(all_matchings) == 105
    tasks += [("all_105_matchings", blocks, (0, 1, 5, 6), None) for blocks in all_matchings]
    by_name = {}
    process_ids = set()
    total_cpu_seconds = 0
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for name, pid, seconds, result in pool.map(task, tasks):
            by_name.setdefault(name, Counter()).update(result)
            process_ids.add(pid)
            total_cpu_seconds += seconds
    examples = explicit_examples()
    totals = Counter()
    for value in by_name.values():
        totals.update(value)
    result = {
        "platform": platform.platform(), "python": platform.python_version(),
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "arithmetic": "Integer Walsh, integer GF(2) ranks and exact rational weights.",
        "workers_configured": args.workers, "distinct_workers_observed": len(process_ids),
        "worker_process_cpu_seconds": total_cpu_seconds,
        "cases_by_geometries": by_name, "totals": totals, "examples": examples,
        "elapsed_seconds": time.perf_counter() - start,
    }
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(totals, sort_keys=True))


if __name__ == "__main__":
    main()
