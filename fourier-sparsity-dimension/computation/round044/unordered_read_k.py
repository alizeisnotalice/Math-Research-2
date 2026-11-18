#!/usr/bin/env python3
"""Independent isolation checks for read-K trees without repeated labels along a path."""

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import hashlib
from itertools import product
import json
from pathlib import Path
import platform
import time

from read_k_layers import L, LARGE, SMALL, Y, Tree, fwht, node, spectra


def compact(labelled, bits):
    return sum(((labelled >> bit) & 1) << i for i, bit in enumerate(bits))


def plan(tree):
    multiplicity = Counter()
    groups = {}
    for v, (j, _, _, path) in enumerate(tree.nodes):
        if j >= 0:
            assert len({label for label, _ in path} | {j}) == len(path) + 1
            groups.setdefault(j, []).append(v)
            multiplicity.update(a for a, _ in tree.fibres[v])
    assert max(multiplicity.values()) <= (1 << (tree.K - 1))
    plans = {}
    for j, same in groups.items():
        j_bit = 1 << j
        for v in same:
            own_mask, own_value = tree.path_conditions[v]
            fixed = 0
            for w in same:
                if w == v:
                    continue
                other_mask, other_value = tree.path_conditions[w]
                lhs, rhs = next(
                    (lhs, rhs) for lhs, rhs in zip(tree.nodes[v][3], tree.nodes[w][3])
                    if lhs != rhs
                )
                assert lhs[0] == rhs[0] and lhs[1] != rhs[1]
                blocker = 1 << lhs[0]
                assert blocker & own_mask & other_mask & (own_value ^ other_value)
                fixed |= blocker
            c = fixed.bit_count()
            assert c <= len(same) - 1 <= tree.K - 1
            free = [bit for bit in range(tree.label_count) if (own_mask & ~fixed) >> bit & 1]
            used = own_mask | j_bit
            physical = [
                (x, compact(x, free), not (x ^ own_value) & own_mask)
                for x in range(tree.selector_size)
                if not (x & j_bit) and not (x ^ own_value) & fixed
            ]
            freq = [
                (a, compact(a, free), 1 - 2 * ((a & fixed & own_value).bit_count() & 1))
                for a, _ in tree.fibres[v]
            ]
            plans[v] = j_bit, c, used, free, physical, freq
    return plans, multiplicity


def audit(tree, plans, family, counts, witnesses):
    w, root, differences, profiles, marks = spectra(tree, family)
    haar = root + [0] * (len(w) - Y)
    signs = [[1 - 2 * ((mask >> y) & 1) for y in range(Y)] for mask in family]
    for v, profile in profiles.items():
        for a, sign in tree.fibres[v]:
            for b, val in enumerate(profile):
                size = len(tree.fibres[v])
                assert val % size == 0
                haar[a * Y + b] += sign * (val // size)
    assert haar == w
    counts["full_Walsh_coefficients_checked"] += len(w)
    for v, profile in profiles.items():
        j_bit, c, used, free, physical, freq = plans[v]
        band_size = sum(w[a * Y + b] != 0 for a, _, _ in freq for b in range(Y))
        isolated_raw = [0] * (Y << len(free))
        checked_cell = {}
        for x, x_sub, own in physical:
            plus, minus = signs[tree.routes[x]], signs[tree.routes[x ^ j_bit]]
            for b, val in enumerate(a - b for a, b in zip(plus, minus)):
                isolated_raw[x_sub * Y + b] += val
            if x_sub in checked_cell:
                assert checked_cell[x_sub] == own
            checked_cell[x_sub] = own
        outside = 1 << (tree.label_count - used.bit_count())
        for x_sub, own in checked_cell.items():
            for b, val in enumerate(differences[v]):
                assert tree.selector_size * isolated_raw[x_sub * Y + b] == (
                    2 * outside * val if own else 0
                )
        isolated = fwht(isolated_raw)
        restricted_from_full_spectrum = [0] * len(isolated)
        for a, a_sub, phase in freq:
            for b in range(Y):
                restricted_from_full_spectrum[a_sub * Y + b] += phase * w[a * Y + b]
        assert restricted_from_full_spectrum == [val << c for val in isolated]
        actual_iso = sum(val != 0 for val in isolated)
        predicted_iso = (1 << len(free)) * sum(val != 0 for val in profile)
        assert actual_iso == predicted_iso and band_size >= actual_iso
        assert band_size * (1 << (tree.K - 1)) >= (
            len(tree.fibres[v]) * sum(val != 0 for val in profile)
        )
        counts["marked_nodes_isolated"] += 1
        counts["tight_global_K_band_bounds"] += band_size * (1 << (tree.K - 1)) == (
            len(tree.fibres[v]) * sum(val != 0 for val in profile)
        )
        same_marked = len(marks[j_bit.bit_length() - 1])
        if band_size * (1 << (same_marked - 1)) < (
            len(tree.fibres[v]) * sum(val != 0 for val in profile)
        ):
            counts["marked_only_comparison_failures"] += 1
            if "marked_only_comparison_failure" not in witnesses:
                witnesses["marked_only_comparison_failure"] = {
                    "leaf_truth_masks": family,
                    "node_index": v,
                    "H_v_actual_sparsity": band_size,
                    "depth": len(tree.fibres[v]).bit_length() - 1,
                    "same_label_marked_nodes": same_marked,
                }


def audit_template(args):
    name, shape, library = args
    tree = Tree(shape, ordered=False)
    plans, multiplicity = plan(tree)
    counts = Counter()
    witnesses = {}
    for masks in product(library, repeat=tree.leaf_count):
        audit(tree, plans, masks, counts, witnesses)
        counts["families"] += 1
    return name, {
        "shape": shape,
        "physical_read_K": tree.K,
        "leaf_library": library,
        "selector_band_multiplicities": dict(multiplicity),
        "maximum_overlap": max(multiplicity.values()),
        "counts": dict(counts),
        "witnesses": witnesses,
    }


def not_just_marked_example():
    shape = node(0, node(2, node(1), node(1)), node(1, node(2), node(2)))
    masks = (0, 15, 15, 0, 0, 15, 15, 0)
    tree = Tree(shape, ordered=False)
    plans, multiplicity = plan(tree)
    counts = Counter()
    findings = {}
    audit(tree, plans, masks, counts, findings)
    actual, _, _, profiles, marked_groups = spectra(tree, masks)
    v = next(i for i, (j, _, _, p) in enumerate(tree.nodes) if j == 1 and p == ((0, 0), (2, 0)))
    dormant = next(i for i, (j, _, _, p) in enumerate(tree.nodes) if j == 1 and p == ((0, 1),))
    assert v in profiles and dormant not in profiles
    assert (tree.routes[1], tree.routes[3]) == (4, 6)
    assert masks[tree.routes[1]] == 0 and masks[tree.routes[3]] == 15
    band_size = sum(actual[a * Y + b] != 0 for a, _ in tree.fibres[v] for b in range(Y))
    marked_same_label = len(marked_groups[1])
    depth = len(tree.nodes[v][3])
    delta_sparsity = sum(coefficient != 0 for coefficient in profiles[v])
    assert (band_size, marked_same_label, depth, delta_sparsity) == (1, 2, 2, 1)
    assert band_size * (1 << (marked_same_label - 1)) < (1 << depth) * delta_sparsity
    assert [a for a, val in enumerate(actual) if val] == [24] and multiplicity[6] == 4
    return {
        "shape": shape,
        "leaf_truth_masks": masks,
        "explicit_Walsh_support": [a for a, val in enumerate(actual) if val],
        "target_node_index": v,
        "target_depth": depth,
        "target_H_v_actual_sparsity": band_size,
        "target_delta_sparsity": delta_sparsity,
        "target_same_label_marked_nodes": marked_same_label,
        "target_same_label_physical_nodes": sum(j == 1 for j, _, _, _ in tree.nodes),
        "dormant_node_index": dormant,
        "dormant_becomes_active_after_fixing_label_2_to_zero": True,
        "overlap_at_selector_frequency_0b110": multiplicity[6],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--workers", type=int, default=6)
    args = parser.parse_args()
    started = time.perf_counter()
    tasks = [
        ("crossed_orders", node(0, node(2, node(1), node(1)), node(1, node(2), node(2))), LARGE),
        ("ascending_control", node(0, node(1, node(2), node(2)), node(1, node(2), node(2))), LARGE),
        ("descending_shared", node(3, node(0, node(2), L), node(2, L, node(0))), LARGE),
        ("nonuniform_depth", node(2, node(3, node(1), L), node(1, node(0), L)), LARGE),
        ("reversal_with_two_reuses", node(2, node(3, node(0), L), node(0, node(3), L)), LARGE),
        ("two_nested_reversals", node(2, node(1, node(0), L), node(0, node(3), node(3))), LARGE),
        ("simple_reversed_order", node(1, node(0), node(0)), SMALL),
    ]
    ans = {}
    totals = Counter()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for name, result in pool.map(audit_template, tasks):
            ans[name] = result
            totals.update(result["counts"])
            print(name, result["counts"]["families"], flush=True)
    record = {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "parallel_workers": args.workers,
        "arithmetic": "Exact integers; complete finite enumeration.",
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "ordered_helper_sha256": hashlib.sha256(Path(__file__).with_name("read_k_layers.py").read_bytes()).hexdigest(),
        "templates": ans,
        "totals": dict(totals),
        "counterexample_to_marked_only_isolation": not_just_marked_example(),
        "elapsed_seconds": time.perf_counter() - started,
    }
    args.output.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
