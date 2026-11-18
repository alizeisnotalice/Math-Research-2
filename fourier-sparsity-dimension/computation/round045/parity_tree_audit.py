#!/usr/bin/env python3
"""Exact local isolation and Kraft costing for truly dependent parity trees."""

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import hashlib
from itertools import product
import json
from pathlib import Path
import platform
import random
import time


Y = 4
LIBRARY = (0, 15, 3, 1, 2)
SEED = 202611180945
LEAF = None


def node(q, left=LEAF, right=LEAF, offset=0):
    return q, offset, left, right


def parity(x):
    return x.bit_count() & 1


def fwht(values):
    data = list(values)
    step = 1
    while step < len(data):
        for start in range(0, len(data), step * 2):
            for i in range(start, start + step):
                j = i + step
                a, b = data[i], data[j]
                data[i], data[j] = a + b, a - b
        step *= 2
    return data


def span(generators):
    out = {0}
    for v in generators:
        out |= {v ^ w for w in out}
    return out


def solve_gf2(coefficients, rhs, n):
    rows = [q | (bit << n) for q, bit in zip(coefficients, rhs)]
    leading = []
    for column in range(n):
        r = len(leading)
        pivot = next((i for i in range(r, len(rows)) if (rows[i] >> column) & 1), None)
        if pivot is None:
            continue
        rows[r], rows[pivot] = rows[pivot], rows[r]
        for i in range(len(rows)):
            if i != r and (rows[i] >> column) & 1:
                rows[i] ^= rows[r]
        leading.append(column)
    if any(r == 1 << n for r in rows):
        return [], len(leading)
    origin = sum(((row >> n) & 1) << pivot for row, pivot in zip(rows, leading))
    kernel = [
        (1 << bit) | sum(((row >> bit) & 1) << pivot for row, pivot in zip(rows, leading))
        for bit in range(n)
        if bit not in leading
    ]
    return sorted(origin ^ x for x in span(kernel)), len(leading)


class Tree:
    def __init__(self, shape):
        self.nodes = []
        self.leaf_count = 0
        self.compile(shape, ())
        self.internal = [v for v, (q, _, _, _, _) in enumerate(self.nodes) if q]
        self.variables = max((q.bit_length() for q, _, _, _, _ in self.nodes), default=0)
        self.domain = 1 << self.variables
        self.routes = [self.leaf(x) for x in range(self.domain)]
        self.words = {self.nodes[v][0] for v in self.internal}

    def compile(self, shape, path):
        index = len(self.nodes)
        self.nodes.append(None)
        if shape is None:
            self.nodes[index] = 0, 0, self.leaf_count, -1, path
            self.leaf_count += 1
        else:
            q, offset, left, right = shape
            lhs = self.compile(left, path + ((index, q, offset, 0),))
            rhs = self.compile(right, path + ((index, q, offset, 1),))
            self.nodes[index] = q, offset, lhs, rhs, path
        return index

    def leaf(self, x):
        v = 0
        while True:
            q, offset, lhs, rhs, _ = self.nodes[v]
            if not q:
                return lhs
            v = rhs if parity(q & x) ^ offset else lhs

    def blockers(self, v, delta):
        _, _, _, _, path = self.nodes[v]
        blockers = set()
        for u in self.internal:
            q, _, _, _, path_u = self.nodes[u]
            if not parity(q & delta) or u == v or any(i == v for i, *_ in path_u):
                continue
            assert all(i != u for i, *_ in path)
            index = next(i for i, (a, b) in enumerate(zip(path, path_u)) if a != b)
            assert path[index][:3] == path_u[index][:3]
            assert path[index][3] != path_u[index][3]
            blockers.add(index)
        return blockers


def geometry(tree):
    plans = {}
    for v in tree.internal:
        q, offset, _, _, path = tree.nodes[v]
        anc = [form for _, form, _, _ in path]
        depth = len(anc)
        candidates, path_rank = solve_gf2(anc + [q], [0] * depth + [1], tree.variables)
        assert path_rank == depth + 1 and candidates
        assert candidates == [
            x for x in range(tree.domain) if parity(x & q) and all(not parity(x & p) for p in anc)
        ]
        blockers = [(len(cells), delta, cells) for delta in candidates
                    for cells in [tree.blockers(v, delta)]]
        ell, delta, fixed = min(blockers, key=lambda item: (item[0], item[1]))
        free = [i for i in range(depth) if i not in fixed]

        fibres = []
        for mask in range(1 << depth):
            a, fixed_phase, haar_phase, a_sub = q, offset, offset, 0
            for i, (_, form, const, side) in enumerate(path):
                if (mask >> i) & 1:
                    a ^= form
                    haar_phase ^= const ^ side
                    fixed_phase ^= const ^ (side if i in fixed else 0)
            for sub, i in enumerate(free):
                a_sub |= ((mask >> i) & 1) << sub
            fibres.append((a, a_sub, 1 - 2 * fixed_phase, 1 - 2 * haar_phase))
        assert len({a for a, *_ in fibres}) == (1 << depth)
        assert q not in span(anc)

        physical = []
        for x in range(tree.domain):
            bits = [parity(x & p) ^ const for _, p, const, _ in path]
            if parity(x & q) ^ offset or any(bits[i] != path[i][3] for i in fixed):
                continue
            x_sub = sum(bits[i] << sub for sub, i in enumerate(free))
            own = all(bit == side for bit, (_, _, _, side) in zip(bits, path))
            if not own:
                assert tree.routes[x] == tree.routes[x ^ delta]
            physical.append((tree.routes[x], tree.routes[x ^ delta], x_sub, own))
        outside = 1 << (tree.variables - depth - 1)
        assert len(physical) == outside * (1 << len(free))
        plans[v] = {
            "q": q, "depth": depth, "ell": ell, "delta": delta,
            "free_count": len(free), "outside": outside, "physical": physical, "fibres": fibres,
            "lca_positions": sorted(fixed), "path_rank": path_rank, "delta_count": len(candidates),
            "max_ell_over_delta_choices": max(n for n, _, _ in blockers),
        }
    max_ell = max((p["ell"] for p in plans.values()), default=0)
    scale = 1 << max_ell
    weights = Counter()
    multiplicity = Counter()
    for p in plans.values():
        for a, *_ in p["fibres"]:
            weights[a] += 1 << (max_ell - p["ell"])
            multiplicity[a] += 1
    assert all(a and w <= scale for a, w in weights.items())
    return plans, max_ell, weights, multiplicity


def recursive_profiles(tree, signs):
    values = [[] for _ in tree.nodes]
    marks = {}
    for v in reversed(range(len(tree.nodes))):
        q, _, lhs, rhs, _ = tree.nodes[v]
        if not q:
            values[v] = [val * tree.domain for val in signs[lhs]]
            continue
        pairs = list(zip(values[lhs], values[rhs]))
        assert all((a - b) % 2 == 0 for a, b in pairs)
        values[v] = [(a + b) // 2 for a, b in pairs]
        delta = [(a - b) // 2 for a, b in pairs]
        if any(delta):
            marks[v] = delta, fwht(delta)
    return fwht(values[0]), marks


def audit_case(tree, plans, max_ell, family, counts, witnesses):
    signs = [[1 - 2 * ((mask >> y) & 1) for y in range(Y)] for mask in family]
    direct = fwht(value for terminal in tree.routes for value in signs[terminal])
    root, marks = recursive_profiles(tree, signs)
    assembled = root + [0] * (len(direct) - Y)
    cost_num = 0
    scale = 1 << (2 * max_ell)
    for v, (delta_values, delta_spectrum) in marks.items():
        p = plans[v]
        ell, depth, outside, free = p["ell"], p["depth"], p["outside"], p["free_count"]
        nterms = 1 << depth
        physical_odd = [0] * (Y << free)
        projected_band = [0] * (Y << free)
        band_sparsity = 0
        for a, a_sub, restriction_sign, haar_sign in p["fibres"]:
            for b, d in enumerate(delta_spectrum):
                assert d % nterms == 0
                assembled[a * Y + b] += haar_sign * (d // nterms)
                true_value = direct[a * Y + b]
                band_sparsity += true_value != 0
                projected_band[a_sub * Y + b] += restriction_sign * true_value
        own_states = set()
        for plus, minus, x_sub, own in p["physical"]:
            if own:
                own_states.add(x_sub)
            for b, (u, w) in enumerate(zip(signs[plus], signs[minus])):
                physical_odd[x_sub * Y + b] += u - w
        assert len(own_states) == 1
        for x_sub in range(1 << free):
            for b, val in enumerate(delta_values):
                assert physical_odd[x_sub * Y + b] * tree.domain == (
                    2 * outside * val if x_sub in own_states else 0
                )
        isolated = fwht(physical_odd)
        assert projected_band == [val << ell for val in isolated]
        counts["isolated_coefficients_checked"] += len(isolated)
        weight = sum(val != 0 for val in isolated)
        delta_sparsity = sum(val != 0 for val in delta_spectrum)
        assert weight == (1 << (depth - ell)) * delta_sparsity
        assert band_sparsity >= weight
        cost_num += (1 << (depth + 2 * (max_ell - ell))) * delta_sparsity
        counts["marked_nodes_checked"] += 1
        counts["tight_local_isolations"] += band_sparsity == weight
        counts["marked_nodes_with_no_outside_sensitive_query"] += ell == 0
        if ell and band_sparsity < nterms * delta_sparsity and "uncorrected_local_bound_failure" not in witnesses:
            witnesses["uncorrected_local_bound_failure"] = {
                "leaf_masks": family, "node_index": v, "depth": depth, "ell": ell,
                "band_sparsity": band_sparsity, "delta_sparsity": delta_sparsity,
            }
    assert assembled == direct
    counts["full_walsh_coefficients_checked"] += len(direct)
    nonzero_selector = sum(w != 0 for w in direct[Y:])
    assert cost_num <= scale * nonzero_selector
    counts["checks_with_nonzero_marks"] += bool(marks)
    counts["checks_with_no_marks"] += not marks
    counts["tight_nonempty_global_bounds"] += bool(marks) and cost_num == scale * nonzero_selector
    if marks and cost_num == scale * nonzero_selector and "tight_global_bound" not in witnesses:
        witnesses["tight_global_bound"] = {
            "leaf_masks": family, "weighted_cost_numerator": cost_num, "weighted_cost_denominator": scale,
            "nonzero_selector_sparsity": nonzero_selector, "root_data_sparsity": sum(w != 0 for w in root),
        }


def component_check(tree, plans, multiplicity, groups):
    if not groups:
        return {}
    labels = [q for group in groups for q in group]
    assert set(labels) == tree.words and len(labels) == len(set(labels))
    ranks = [len(span(group)).bit_length() - 1 for group in groups]
    assert sum(ranks) == len(span(labels)).bit_length() - 1
    nullities = [len(g) - r for g, r in zip(groups, ranks)]
    component_of = {q: i for i, group in enumerate(groups) for q in group}
    path_nullities = [
        sum(nullities[i] for i in {component_of[q] for _, q, _, _ in path})
        for q, _, _, _, path in tree.nodes if q == 0
    ]
    K = max(Counter(tree.nodes[v][0] for v in tree.internal).values())
    assert max(multiplicity.values()) <= (1 << (max(path_nullities) + K - 1))
    assert all(p["ell"] <= K * (1 + nullities[component_of[p["q"]]]) - 1 for p in plans.values())
    return {"query_groups": groups, "ranks": ranks, "nullities": nullities, "maximum_path_nullity": max(path_nullities)}


def run_template(task):
    name, shape, library, sample_size, groups = task
    start = time.perf_counter()
    tree = Tree(shape)
    plans, max_ell, weights, mult = geometry(tree)
    stats = Counter()
    witnesses = {}
    if sample_size:
        rng = random.Random(SEED)
        sampled = {(0,) * tree.leaf_count, (15,) * tree.leaf_count,
                   tuple(LIBRARY[i % len(LIBRARY)] for i in range(tree.leaf_count))}
        while len(sampled) < sample_size:
            sampled.add(tuple(rng.choice(library) for _ in range(tree.leaf_count)))
        families = sorted(sampled)
    else:
        families = product(library, repeat=tree.leaf_count)
    for family in families:
        audit_case(tree, plans, max_ell, family, stats, witnesses)
        stats["families"] += 1
    return name, {
        "shape": shape,
        "leaf_library": library,
        "sampling_seed": SEED if sample_size else None,
        "physical_selector_variables": tree.variables,
        "distinct_query_directions": sorted(tree.words),
        "query_rank": len(span(tree.words)).bit_length() - 1,
        "maximum_isolation_cost": max_ell,
        "maximum_actual_band_multiplicity": max(mult.values(), default=0),
        "tight_Kraft_frequencies": sorted(a for a, val in weights.items() if val == (1 << max_ell)),
        "component_check": component_check(tree, plans, mult, groups),
        "node_certificates": {
            v: {k: p[k] for k in ("q", "depth", "ell", "delta", "lca_positions", "path_rank",
                                  "delta_count", "max_ell_over_delta_choices")}
            for v, p in plans.items()
        },
        "counts": dict(stats),
        "witnesses": witnesses,
        "seconds": time.perf_counter() - start,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--workers", type=int, default=6)
    args = parser.parse_args()
    start = time.perf_counter()
    a, b, c = node(4, node(8), node(12)), node(16, node(32), node(48)), node(64, node(128), node(192))
    tasks = [
        ("dependent_triad_all_leaves", node(1, node(2), node(3)), tuple(range(16)), 0, [(1, 2, 3)]),
        ("offset_dependent_triad", node(1, node(2, offset=1), node(3), offset=1), LIBRARY, 0, [(1, 2, 3)]),
        ("multiple_dependencies", node(1, node(2, node(4), LEAF), node(3, LEAF, node(5))), LIBRARY, 0, []),
        ("independent_queries_empty_blockers", node(1, node(2), node(4)), LIBRARY, 0, [(1,), (2,), (4,)]),
        ("seven_distinct_dependent_queries", node(1, node(2, node(4), node(5)), node(3, node(6), node(7))), (0, 15, 1), 0, []),
        ("dependent_and_repeated_query", node(1, node(2, node(4), LEAF), node(3, node(4), LEAF)), (0, 15, 1), 0, []),
        ("three_independent_dependent_clusters", node(1, a, node(2, b, c)), LIBRARY, 128,
         [(1,), (2,), (4, 8, 12), (16, 32, 48), (64, 128, 192)]),
        ("leaf_only", LEAF, tuple(range(16)), 0, []),
    ]
    totals = Counter()
    reports = {}
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for name, result in pool.map(run_template, tasks):
            reports[name] = result
            totals.update(result["counts"])
            print(name, result["counts"]["families"], flush=True)
    assert totals["marked_nodes_with_no_outside_sensitive_query"] and totals["checks_with_no_marks"]
    bad_delta, bad_rank = solve_gf2([1, 2, 3], [0, 0, 1], 2)
    assert not bad_delta and bad_rank == 2
    result = {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "arithmetic": "Exact integers; explicit GF(2) elimination and Walsh transforms.",
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "workers": args.workers,
        "totals": dict(totals),
        "templates": reports,
        "excluded_path_dependent_control": {
            "path_query_vectors": [1, 2, 3], "path_rank": bad_rank, "delta_solutions": bad_delta,
        },
        "seconds": time.perf_counter() - start,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
