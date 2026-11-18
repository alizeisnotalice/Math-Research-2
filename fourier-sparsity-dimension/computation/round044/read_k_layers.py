#!/usr/bin/env python3
"""Parallel exact Walsh and isolation checks for strictly ordered read-K trees."""

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import hashlib
from itertools import product
import json
from pathlib import Path
import platform
import time


Y = 4
SMALL = (0, 15, 3, 1, 2)
LARGE = (0, 15, 1)
L = None


def node(label, lhs=L, rhs=L):
    return (label, lhs, rhs)


def rank(points):
    basis = {}
    for p in points:
        while p:
            i = p.bit_length() - 1
            if i not in basis:
                basis[i] = p
                break
            p ^= basis[i]
    return len(basis)


def fwht(values):
    out = list(values)
    d = 1
    while d < len(out):
        for start in range(0, len(out), 2 * d):
            for i in range(start, start + d):
                j = i + d
                a, b = out[i], out[j]
                out[i], out[j] = a + b, a - b
        d *= 2
    return out


class Tree:
    def __init__(self, shape, *, ordered=True):
        self.nodes = []
        self.leaf_count = 0
        self.compile(shape, -1, (), ordered)
        self.label_count = 1 + max(label for label, _, _, _ in self.nodes if label >= 0)
        self.selector_size = 1 << self.label_count
        self.routes = [self.leaf(x) for x in range(self.selector_size)]
        self.K = max(Counter(j for j, _, _, _ in self.nodes if j >= 0).values())
        self.path_conditions = []
        self.fibres = []
        for label, _, _, path in self.nodes:
            bits = sum(1 << j for j, _ in path)
            value = sum(b << j for j, b in path)
            self.path_conditions.append((bits, value))
            if label < 0:
                self.fibres.append([])
                continue
            frequencies = [(1 << label, 1)]
            for j, side in path:
                frequencies += [(f ^ (1 << j), s * (1 - 2 * side)) for f, s in frequencies]
            self.fibres.append(frequencies)

    def compile(self, shape, lower, path, ordered):
        i = len(self.nodes)
        self.nodes.append(None)
        if shape is None:
            self.nodes[i] = (-1, self.leaf_count, -1, path)
            self.leaf_count += 1
            return i
        label, left, right = shape
        if ordered:
            assert label > lower
        lhs = self.compile(left, label, path + ((label, 0),), ordered)
        rhs = self.compile(right, label, path + ((label, 1),), ordered)
        self.nodes[i] = (label, lhs, rhs, path)
        return i

    def leaf(self, x):
        v = 0
        while True:
            j, lhs, rhs, _ = self.nodes[v]
            if j < 0:
                return lhs
            v = rhs if (x >> j) & 1 else lhs


def spectra(tree, masks):
    signs = [[1 - 2 * ((m >> y) & 1) for y in range(Y)] for m in masks]
    truth = [value for leaf in tree.routes for value in signs[leaf]]
    direct = fwht(truth)
    means = [[] for _ in tree.nodes]
    differences = {}
    profiles = {}
    marks = {}
    for v in reversed(range(len(tree.nodes))):
        label, lhs, rhs, _ = tree.nodes[v]
        if label < 0:
            means[v] = [value * tree.selector_size for value in signs[lhs]]
            continue
        pairs = list(zip(means[lhs], means[rhs]))
        assert all((a - b) % 2 == 0 for a, b in pairs)
        means[v] = [(a + b) // 2 for a, b in pairs]
        delta = [(a - b) // 2 for a, b in pairs]
        if any(delta):
            differences[v] = delta
            profiles[v] = fwht(delta)
            marks.setdefault(label, []).append(v)
    return direct, fwht(means[0]), differences, profiles, marks


def audit(tree, masks, counts, witnesses):
    direct, root_profile, differences, profiles, marks = spectra(tree, masks)
    reconstructed = root_profile + [0] * (len(direct) - Y)
    active_formal = set()
    formal_support = sum(w != 0 for w in root_profile)
    for j, group in marks.items():
        for v in group:
            profile = profiles[v]
            fibres = tree.fibres[v]
            formal_support += sum(w != 0 for w in profile) * len(fibres)
            assert all(w % len(fibres) == 0 for w in profile)
            for a, sign in fibres:
                assert a.bit_length() == j + 1
                for b, w in enumerate(profile):
                    reconstructed[a * Y + b] += sign * (w // len(fibres))
            active_formal.add(v)
    assert direct == reconstructed
    counts["compared_coefficients"] += len(direct)

    def witness(key, **extra):
        if key in witnesses:
            return
        witnesses[key] = {
            "leaf_truth_masks": [hex(m) for m in masks],
            "Walsh_support": [i for i, w in enumerate(direct) if w],
            **extra,
        }

    s = sum(w != 0 for w in direct)
    assert formal_support <= tree.K * (1 << (tree.K - 1)) * s
    selector_support = [xi // Y for xi, w in enumerate(direct) if w]
    support_of_selector_bits = 0
    for xi in selector_support:
        support_of_selector_bits |= xi
    if rank(selector_support) < support_of_selector_bits.bit_count():
        counts["projection_rank_less_than_effective_coordinate_count"] += 1
        witness("projection_rank_loss", rank=rank(selector_support),
                active_bits=support_of_selector_bits.bit_count())
    counts["same_label_layers"] += sum(len(group) > 1 for group in marks.values())
    for j, group in sorted(marks.items()):
        lo = (1 << j) * Y
        hi = 2 * lo
        sj = sum(w != 0 for w in direct[lo:hi])
        naive_layer_size = sum(
            len(tree.fibres[v]) * sum(w != 0 for w in profiles[v]) for v in group
        )
        if sj < naive_layer_size:
            counts["layer_collision_losses"] += 1
        m = len(group)
        for v in group:
            own_mask, own_value = tree.path_conditions[v]
            depth = len(tree.nodes[v][3])
            fixed_mask = 0
            for w in group:
                if w == v:
                    continue
                other_mask, other_value = tree.path_conditions[w]
                conflict = own_mask & other_mask & (own_value ^ other_value)
                assert conflict
                fixed_mask |= conflict & -conflict
            c = fixed_mask.bit_count()
            assert c <= min(depth, m - 1)
            restricted = []
            for x in range(1 << j):
                if (x ^ own_value) & fixed_mask:
                    continue
                participants = [
                    w for w in group if not ((x ^ tree.path_conditions[w][1]) & tree.path_conditions[w][0])
                ]
                assert participants == ([v] if (x ^ own_value) & own_mask == 0 else [])
                restricted.extend(differences[v] if participants else [0] * Y)
            isolated_support = sum(w != 0 for w in fwht(restricted))
            term_support = sum(w != 0 for w in profiles[v])
            assert isolated_support == (1 << (depth - c)) * term_support
            assert sj >= isolated_support
            assert sj * (1 << (m - 1)) >= (1 << depth) * term_support
            counts["marked_isolations_checked"] += 1
            counts["tight_isolation_bounds"] += sj == isolated_support
            if m > 1:
                tight = sj * (1 << (m - 1)) == (1 << depth) * term_support
                counts["tight_read_K_bounds_with_repetition"] += tight
                if tight:
                    witness("sharp_repeated_layer", j=j, marks=m, depth=depth, blockers=c,
                            layer_sparsity=sj, delta_sparsity=term_support)
                if sj < (1 << depth) * term_support:
                    counts["individual_uncorrected_bound_failures"] += 1
                    witness("uncorrected_bound_failure", j=j, marks=m, depth=depth, blockers=c,
                            layer_sparsity=sj, delta_sparsity=term_support)


def one_template(task):
    name, shape, library = task
    tree = Tree(shape)
    counts = Counter()
    witnesses = {}
    for masks in product(library, repeat=tree.leaf_count):
        audit(tree, masks, counts, witnesses)
        counts["families"] += 1
    return name, {
        "shape": shape,
        "leaf_library": [hex(m) for m in library],
        "label_count": tree.label_count,
        "physical_read_K": tree.K,
        "counts": dict(counts),
        "witnesses": witnesses,
    }


def negative_controls():
    examples = [
        ("unordered_suppresses_entire_lower_layer", node(1, node(0), node(0)), (0, 15, 15, 0), 0),
        ("repeating_label_on_path_marks_unreachable_edge", node(0, node(0), L), (0, 15, 0), 0),
    ]
    result = {}
    for name, shape, masks, target_label in examples:
        tree = Tree(shape, ordered=False)
        w, _, _, profiles, marks = spectra(tree, masks)
        lo = Y << target_label
        assert all(entry == 0 for entry in w[lo:2 * lo])
        assert marks[target_label]
        result[name] = {
            "shape": shape,
            "leaf_truth_masks": masks,
            "Walsh_support": [i for i, c in enumerate(w) if c],
            "empty_layer": target_label,
            "formal_marked_profiles": [profiles[v] for v in marks[target_label]],
        }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    started = time.perf_counter()
    tasks = [
        ("read_two_balanced", node(0, node(1), node(1)), SMALL),
        ("read_two_depth_mismatch", node(0, node(3), node(1, node(3), L)), SMALL),
        ("read_three_comb", node(0, node(3), node(1, node(3), node(3))), LARGE),
        ("read_four_comb", node(0, node(3), node(1, node(3), node(2, node(3), node(3)))), LARGE),
        ("read_four_balanced", node(0, node(1, node(2), node(2)), node(1, node(2), node(2))), LARGE),
        ("two_distinct_colliding_layers", node(0, node(1, node(3), L), node(1, L, node(3))), LARGE),
        ("read_one_control", node(0, node(1), node(2)), SMALL),
        ("sparse_labels_and_mixed_depths", node(0, node(4), node(2, node(3), node(4))), LARGE),
    ]
    output = {}
    totals = Counter()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for name, result in pool.map(one_template, tasks):
            output[name] = result
            totals.update(result["counts"])
            print(name, result["counts"]["families"], flush=True)
    assert totals["tight_read_K_bounds_with_repetition"]
    assert totals["individual_uncorrected_bound_failures"]
    record = {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "parallel_workers": args.workers,
        "arithmetic": "Exact integers and complete finite enumeration.",
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "templates": output,
        "totals": dict(totals),
        "negative_controls": negative_controls(),
        "elapsed_seconds": time.perf_counter() - started,
    }
    args.output.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
