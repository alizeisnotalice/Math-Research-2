#!/usr/bin/env python3
"""Exact 2x2-block / 4x2-gluing audit for the peak-budget lemma."""

from collections import Counter
from itertools import combinations, permutations
import json
from pathlib import Path


def graph(width, height):
    adj = []
    for y in range(height):
        for x in range(width):
            neighbours = []
            for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                xx, yy = x + dx, y + dy
                if 0 <= xx < width and 0 <= yy < height:
                    neighbours.append(yy * width + xx)
            adj.append(tuple(neighbours))
    return adj


def peaks(order, adj):
    position = [0] * len(order)
    for i, v in enumerate(order):
        position[v] = i  # high-to-low order
    return tuple(v for v, ns in enumerate(adj) if all(position[v] < position[u] for u in ns))


def interleave(left, right):
    """All binary interleavings, retaining each word's internal order."""
    total = len(left) + len(right)
    for left_slots in combinations(range(total), len(left)):
        left_slots = set(left_slots)
        i = j = 0
        word = []
        for place in range(total):
            if place in left_slots:
                word.append(left[i])
                i += 1
            else:
                word.append(right[j])
                j += 1
        yield tuple(word)


def main():
    small_adj = graph(2, 2)
    words = list(permutations(range(4)))
    by_peak_count = Counter(len(peaks(word, small_adj)) for word in words)
    two_peak_words = [word for word in words if len(peaks(word, small_adj)) == 2]

    # Left block sits at x=0,1; right block at x=2,3 in a 4x2 rectangle.
    full_adj = graph(4, 2)
    global_k = Counter()
    interface_peak_budget_violations = 0
    for left in two_peak_words:
        left_full = tuple((v % 2) + 4 * (v // 2) for v in left)
        for right in two_peak_words:
            right_full = tuple(2 + (v % 2) + 4 * (v // 2) for v in right)
            for word in interleave(left_full, right_full):
                peak_set = peaks(word, full_adj)
                global_k[len(peak_set)] += 1
                # Every local peak away from the interface must remain global.
                left_local = {(v % 2) + 4 * (v // 2) for v in peaks(left, small_adj)}
                right_local = {2 + (v % 2) + 4 * (v // 2) for v in peaks(right, small_adj)}
                non_interface = ((left_local - {1, 5}) | (right_local - {2, 6}))
                if not non_interface.issubset(set(peak_set)):
                    interface_peak_budget_violations += 1

    assert interface_peak_budget_violations == 0

    # A whole K=2 configuration can put both non-interface peaks on one side.
    same_side_example = None
    for order in permutations(range(8)):
        peak_set = peaks(order, full_adj)
        if len(peak_set) == 2 and all((v % 4) <= 1 for v in peak_set):
            same_side_example = {"high_to_low_order": list(order), "peaks": list(peak_set)}
            break
    assert same_side_example is not None
    answer = {
        "two_by_two_peak_histogram": dict(sorted(by_peak_count.items())),
        "two_by_two_A2": len(two_peak_words),
        "two_K2_blocks_all_70_interleavings_global_peak_histogram": dict(sorted(global_k.items())),
        "total_pairs_times_interleavings": len(two_peak_words) ** 2 * 70,
        "two_noninterface_peaks_on_one_side_of_cut_after_column_3": same_side_example,
        "status": "passed",
    }
    Path("peak-budget-audit.json").write_text(json.dumps(answer, indent=2) + "\n")
    print(json.dumps(answer, indent=2))


if __name__ == "__main__":
    main()
