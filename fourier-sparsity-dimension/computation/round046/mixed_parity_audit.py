#!/usr/bin/env python3
"""Exact Walsh audit of mixed shared affine and dependent private parity trees."""

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
import hashlib
from itertools import product
import json
from pathlib import Path
import platform
import runpy
import time


W_BITS = 2
W = 1 << W_BITS
COMMON_LEAVES = (0, 15, 3, 1, 2)
L = None


def configure(path):
    global fwht, solve_gf2, parity, span
    helpers = runpy.run_path(str(path))
    fwht, solve_gf2 = helpers["fwht"], helpers["solve_gf2"]
    parity, span = helpers["parity"], helpers["span"]


def P(q, b, left=L, right=L, const=0):
    return q, b, const, left, right


def S(a, left=L, right=L, const=0):
    return 0, a, const, left, right


class Tree:
    def __init__(self, shape):
        self.nodes = []
        self.leaf_count = 0
        self.compile(shape, ())
        self.privates = [v for v, (q, _, _, _, _, _) in enumerate(self.nodes) if q > 0]
        self.x_bits = max((q.bit_length() for q, _, _, _, _, _ in self.nodes if q > 0), default=0)
        self.X = 1 << self.x_bits
        self.routes = [self.route(x, w) for x in range(self.X) for w in range(W)]

    def compile(self, shape, path):
        i = len(self.nodes)
        self.nodes.append(None)
        if shape is None:
            self.nodes[i] = -1, 0, 0, self.leaf_count, -1, path
            self.leaf_count += 1
        else:
            q, a, b, left, right = shape
            left_index = self.compile(left, path + ((i, q, a, b, 0),))
            right_index = self.compile(right, path + ((i, q, a, b, 1),))
            self.nodes[i] = q, a, b, left_index, right_index, path
        return i

    def route(self, x, w):
        v = 0
        while True:
            q, a, b, left, right, _ = self.nodes[v]
            if q == -1:
                return left
            v = right if (parity(q & x) ^ parity(a & w) ^ b) else left

    def blockers(self, v, delta):
        path = self.nodes[v][-1]
        positions = set()
        for u in self.privates:
            q, _, _, _, _, upath = self.nodes[u]
            if u == v or not parity(q & delta) or any(a[0] == v for a in upath):
                continue
            assert all(a[0] != u for a in path)
            first = next(i for i, (a, b) in enumerate(zip(path, upath)) if a != b)
            assert path[first][:-1] == upath[first][:-1] and path[first][-1] != upath[first][-1]
            positions.add(first)
        return positions


def affine_parameterization(points):
    if not points:
        return [], []
    origin = points[0]
    basis = []
    kernel = {0}
    ordered = [origin]
    for w in points:
        if w ^ origin not in kernel:
            new = w ^ origin
            basis.append(new)
            kernel |= {x ^ new for x in kernel}
            ordered.extend(x ^ new for x in ordered[:])
    assert set(ordered) == set(points)
    return ordered, basis


def shared_space(path, positions):
    equations = [(a, b ^ side) for i, (_, q, a, b, side) in enumerate(path) if not q and i in positions]
    solved, rank = solve_gf2([a for a, _ in equations], [c for _, c in equations], W_BITS)
    assert solved == [w for w in range(W) if all(parity(a & w) == b for a, b in equations)]
    if solved:
        assert len(solved) == (W >> rank)
    return solved, rank


def geometry(tree):
    plans = {}
    for v in tree.privates:
        q, a, b, _, _, path = tree.nodes[v]
        anc_private = [i for i, entry in enumerate(path) if entry[1]]
        k = len(anc_private)
        private_forms = [path[i][1] for i in anc_private]
        solutions, path_rank = solve_gf2(private_forms + [q], [0] * k + [1], tree.x_bits)
        assert path_rank == k + 1 and solutions
        assert solutions == [
            x for x in range(tree.X) if parity(q & x) and not any(parity(form & x) for form in private_forms)
        ]
        options = [(len(pos), delta, pos) for delta in solutions for pos in [tree.blockers(v, delta)]]
        ell, delta, positions = min(options, key=lambda entry: (entry[0], entry[1]))
        private_fixed = [i for i in anc_private if i in positions]
        free = [i for i in anc_private if i not in positions]
        shared_lcas = ell - len(private_fixed)
        reachable, _ = shared_space(path, set(range(len(path))))
        blocked, shared_rank = shared_space(path, positions)
        w_coordinates, w_basis = affine_parameterization(blocked)
        assert set(reachable) <= set(blocked) and shared_rank <= shared_lcas
        bands = []
        for mask in range(1 << k):
            private_freq, data_modulation = q, a
            haar_phase, restriction_phase = b, b
            frequency_free = 0
            for j, i in enumerate(anc_private):
                if (mask >> j) & 1 == 0:
                    continue
                _, normal, shift, const, branch = path[i]
                private_freq ^= normal
                data_modulation ^= shift
                haar_phase ^= const ^ branch
                restriction_phase ^= const ^ (branch if i in private_fixed else 0)
                if i in free:
                    frequency_free |= 1 << free.index(i)
            data_restriction = []
            if w_coordinates:
                for beta in range(W):
                    adjusted = beta ^ data_modulation
                    w_frequency = sum(parity(adjusted & direction) << j for j, direction in enumerate(w_basis))
                    sign = 1 - 2 * (restriction_phase ^ parity(adjusted & w_coordinates[0]))
                    data_restriction.append((w_frequency, sign))
            bands.append((private_freq, data_modulation, 1 - 2 * haar_phase, frequency_free, data_restriction))
        assert len({h for h, *_ in bands}) == (1 << k)
        physical = []
        outside = 1 << (tree.x_bits - k - 1)
        for w_sub, w in enumerate(w_coordinates):
            for x in range(tree.X):
                if parity(q & x) ^ parity(a & w) ^ b:
                    continue
                answers = [parity(px & x) ^ parity(pw & w) ^ pb for _, px, pw, pb, _ in path]
                if any(answers[i] != path[i][-1] for i in private_fixed):
                    continue
                own = all(answer == entry[-1] for answer, entry in zip(answers, path))
                left, right = tree.routes[x * W + w], tree.routes[(x ^ delta) * W + w]
                if not own:
                    assert left == right
                free_state = sum(answers[i] << j for j, i in enumerate(free))
                physical.append((w, w_sub, free_state, left, right, own))
        assert len(physical) == len(blocked) * (1 << len(free)) * outside
        plans[v] = {
            "private_depth": k, "ell": ell, "delta": delta,
            "delta_count": len(solutions), "blocked_ancestor_positions": sorted(positions),
            "path_rank": path_rank, "shared_blocker_rank": shared_rank,
            "shared_lcas": shared_lcas, "private_lcas": len(private_fixed),
            "reachable_w": reachable, "blocked_w": w_coordinates,
            "free_count": len(free), "outside": outside, "physical": physical, "bands": bands,
        }
    maxell = max((p["ell"] for p in plans.values()), default=0)
    scale = 1 << maxell
    charging = Counter()
    for p in plans.values():
        for h, *_ in p["bands"]:
            charging[h] += 1 << (maxell - p["ell"])
    assert all(h and weight <= scale for h, weight in charging.items())
    return plans, maxell, charging


def recursive_haar(tree, plans, signs):
    mean = [[] for _ in tree.nodes]
    marks = {}
    for v in reversed(range(len(tree.nodes))):
        q, a, const, left, right, _ = tree.nodes[v]
        if q == -1:
            mean[v] = [x * tree.X for x in signs[left]]
        elif q == 0:
            mean[v] = [mean[right if parity(a & w) ^ const else left][w] for w in range(W)]
        else:
            pairs = list(zip(mean[left], mean[right]))
            assert all((u - w) % 2 == 0 for u, w in pairs)
            mean[v] = [(u + w) // 2 for u, w in pairs]
            profile = [
                (u - t) // 2 if w in plans[v]["reachable_w"] else 0
                for w, (u, t) in enumerate(pairs)
            ]
            if any(profile):
                marks[v] = profile, fwht(profile)
    return fwht(mean[0]), marks


def audit(tree, plans, maxell, family, counts, witnesses):
    signs = [[1 - 2 * ((mask >> y) & 1) for y in range(W)] for mask in family]
    real = fwht(signs[leaf][i % W] for i, leaf in enumerate(tree.routes))
    root, marks = recursive_haar(tree, plans, signs)
    reconstructed = root + [0] * (len(real) - W)
    global_numerator = 0
    rank_sensitive_numerator = 0
    scale = (1 << (2 * maxell)) * W
    for v, (values, profile) in marks.items():
        p = plans[v]
        reach = len(p["reachable_w"])
        blocked = len(p["blocked_w"])
        assert reach > 0 and blocked > 0
        k, ell, free = p["private_depth"], p["ell"], p["free_count"]
        rank_cost = p["private_lcas"] + p["shared_blocker_rank"]
        size_layer = 1 << k
        iso_size = blocked << free
        raw = [0] * iso_size
        from_full = [0] * iso_size
        own_states = set()
        for w, w_sub, z_free, left, right, own in p["physical"]:
            z = z_free * blocked + w_sub
            raw[z] += signs[left][w] - signs[right][w]
            if own:
                own_states.add(z)
        assert len(own_states) == reach
        for z, val in enumerate(raw):
            w = p["blocked_w"][z % blocked]
            assert val * tree.X == (2 * p["outside"] * values[w] if z in own_states else 0)
        band_sparsity = 0
        for private_freq, modulation, haar_sign, frequency_free, quotient in p["bands"]:
            for beta, (w_sub, sign) in enumerate(quotient):
                val = profile[beta ^ modulation]
                assert val % size_layer == 0
                reconstructed[private_freq * W + beta] += haar_sign * (val // size_layer)
                true_val = real[private_freq * W + beta]
                band_sparsity += true_val != 0
                from_full[frequency_free * blocked + w_sub] += true_val * sign
        isolated = fwht(raw)
        assert from_full == [val << rank_cost for val in isolated]
        s_iso = sum(val != 0 for val in isolated)
        assert s_iso == (1 << free) * sum(bool(x) for x in fwht(values[w] for w in p["blocked_w"]))
        assert band_sparsity >= s_iso
        assert s_iso * reach >= (1 << free) * blocked
        assert (band_sparsity * reach) << rank_cost >= W << k
        assert (band_sparsity * reach) << ell >= W << k
        global_numerator += (W // reach) * W * (1 << (k + 2 * (maxell - ell)))
        rank_sensitive_numerator += (W // reach) * W * (1 << (k + 2 * maxell - rank_cost - ell))
        counts["isolated_coefficients"] += len(isolated)
        counts["marked_private_nodes"] += 1
        counts["marked_nodes_with_shared_lca"] += p["shared_lcas"] > 0
        counts["marked_nodes_with_both_LCA_types"] += p["shared_lcas"] > 0 and p["private_lcas"] > 0
        counts["marked_nodes_with_dependent_shared_blockers"] += p["shared_blocker_rank"] < p["shared_lcas"]
        counts["tight_local_bounds"] += (band_sparsity * reach) << ell == W << k
        if ((band_sparsity * reach) << p["private_lcas"]) < (W << k):
            counts["ignore_shared_LCA_failures"] += 1
            if "drop_shared_LCA" not in witnesses:
                witnesses["drop_shared_LCA"] = {
                    "leaf_truth_masks": family, "node": v, "private_depth": k,
                    "reach_numerator": reach, "reach_denominator": W,
                    "private_lcas": p["private_lcas"], "shared_lcas": p["shared_lcas"],
                    "actual_band_size": band_sparsity,
                }
    assert reconstructed == real
    nonzero_selector = sum(val != 0 for val in real[W:])
    assert global_numerator <= rank_sensitive_numerator <= scale * nonzero_selector
    counts["root_walsh_coefficients"] += len(real)
    counts["empty_mark_cases"] += not marks
    counts["nonempty_tight_global_bounds"] += bool(marks) and global_numerator == scale * nonzero_selector
    if marks and global_numerator == scale * nonzero_selector and "global_tight" not in witnesses:
        witnesses["global_tight"] = {
            "leaf_truth_masks": family, "weighted_cost": str(Fraction(global_numerator, scale)),
            "sparsity": sum(v != 0 for v in real), "nonzero_private_sparsity": nonzero_selector,
        }


def run_template(spec):
    name, shape, library = spec
    tree = Tree(shape)
    plans, maxell, weights = geometry(tree)
    stats = Counter()
    witnesses = {}
    for masks in product(library, repeat=tree.leaf_count):
        audit(tree, plans, maxell, masks, stats, witnesses)
        stats["leaf_families"] += 1
    stats["static_dead_private_nodes"] = sum(not p["reachable_w"] for p in plans.values())
    return name, {
        "shape": shape, "leaf_truth_library": library, "private_variables": tree.x_bits,
        "counts": dict(stats), "witnesses": witnesses,
        "max_ell": maxell, "scaled_Kraft_weights": dict(weights),
        "node_plans": {
            v: {key: val for key, val in p.items() if key not in {"physical", "bands"}}
            for v, p in plans.items()
        }
    }


def nonaffine_counterexample():
    # One shared AND, then two private occurrences which both output the same character.
    table = [1 - 2 * x for x in (0, 1) for _ in range(W)]
    w = fwht(table)
    support = [i for i, a in enumerate(w) if a]
    mass = [Fraction(3, 4), Fraction(1, 4)]
    local = [Fraction(1, 2) / p for p in mass]
    global_cost = sum(Fraction(1, 4) / p for p in mass)
    assert support == [4] and local[-1] == 2 and global_cost > len(support)
    return {
        "test": "Shared AND then the same private character in both branches",
        "Walsh_support": support, "private_blocker_count": 0, "shared_blocker_count": 1,
        "shared_masses": [str(p) for p in mass],
        "unsupported_local_bounds": [str(a) for a in local],
        "unsupported_global_bound": str(global_cost),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("round045_helper", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--workers", type=int, default=6)
    args = parser.parse_args()
    configure(args.round045_helper)
    start = time.perf_counter()
    tasks = [
        ("dependent_repeated_private_under_shared_split",
         S(1, P(1, 2), P(1, 2)), tuple(range(16))),
        ("triad_parities_affine_twists",
         P(1, 1, P(2, 2, const=1), P(3, 2), const=1), tuple(range(16))),
        ("alternating_private_shared_private",
         P(1, 2, S(1, P(2, 2), L), S(1, L, P(3, 2), const=1)), COMMON_LEAVES),
        ("repeated_shared_direction_dead_private",
         S(1, S(1, P(1, 0), P(1, 2)), P(1, 2)), COMMON_LEAVES),
        ("shared_and_private_blockers",
         S(1, P(1, 2, P(2, 1), P(2, 0)), P(2, 2)), COMMON_LEAVES),
        ("shared_only", S(1, S(2, L, L), L), COMMON_LEAVES),
        ("leaf_only", L, tuple(range(16))),
    ]
    totals = Counter()
    result = {}
    with ProcessPoolExecutor(max_workers=args.workers, initializer=configure,
                             initargs=(args.round045_helper,)) as pool:
        for name, values in pool.map(run_template, tasks):
            totals.update(values["counts"])
            result[name] = values
            print(name, values["counts"]["leaf_families"], flush=True)
    assert totals["marked_nodes_with_both_LCA_types"] and totals["marked_nodes_with_dependent_shared_blockers"]
    assert totals["ignore_shared_LCA_failures"] and totals["static_dead_private_nodes"]
    out = {
        "python": platform.python_version(), "platform": platform.platform(),
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "round045_helper_sha256": hashlib.sha256(args.round045_helper.read_bytes()).hexdigest(),
        "arithmetic": "Exact GF(2) matrices and integer Walsh transforms; no tolerance or random seed.",
        "workers": args.workers, "templates": result, "totals": dict(totals),
        "nonaffine_shared_counterexample": nonaffine_counterexample(),
        "elapsed_seconds": time.perf_counter() - start,
    }
    args.output.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
