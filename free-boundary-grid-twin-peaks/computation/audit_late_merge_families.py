#!/usr/bin/env python3
"""Exact finite checks for the parallel-path first-merge configuration."""

from fractions import Fraction
from itertools import permutations
from math import comb, factorial


def neighbors(v, n):
    x, y = v
    for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        q = (x + dx, y + dy)
        if 0 <= q[0] < n and 0 <= q[1] < n:
            yield q


def check(n):
    A = {(0, j) for j in range(n)}
    B = {(2, j) for j in range(n)}
    w = (1, 0)
    assert not any(v in B for u in A for v in neighbors(u, n))
    assert (0, 0) in set(neighbors(w, n))
    assert (2, 0) in set(neighbors(w, n))

    weights = [comb(n - 1, r) for r in range(n)]
    total = sum(weights)
    assert total == 1 << (n - 1)
    mean = sum(Fraction(r * weight, total) for r, weight in enumerate(weights))
    variance = sum(Fraction(weight, total) * (r - mean) ** 2
                   for r, weight in enumerate(weights))
    assert mean == Fraction(n - 1, 2)
    assert variance == Fraction(n - 1, 4)
    # Horizontal root displacement is exactly 2; vertical coordinates are iid.
    squared_distance = 4 + 2 * variance
    assert squared_distance == Fraction(n + 7, 2)
    return {"n": n, "mean_y": str(mean), "var_y": str(variance),
            "mean_separation_squared": 4,
            "conditional_squared_distance": str(squared_distance)}


def brute_path_counts(n):
    """Independent permutation audit of (5), kept to n<=8."""
    counts = [0] * n
    for labels in permutations(range(n)):
        peaks = []
        for r, value in enumerate(labels):
            if (r == 0 or value > labels[r - 1]) and (r == n - 1 or value > labels[r + 1]):
                peaks.append(r)
        if len(peaks) == 1:
            counts[peaks[0]] += 1
    assert counts == [comb(n - 1, r) for r in range(n)], (n, counts)


def concentric_state(n):
    """Exact geometry and boundary-second-moment audit for (10)--(15)."""
    assert n >= 5
    vertices = {(x, y) for x in range(n) for y in range(n)}
    A = {v for v in vertices if min(v[0], v[1], n - 1 - v[0], n - 1 - v[1]) == 0}
    B = {(x, y) for x in range(2, n - 2) for y in range(2, n - 2)}
    assert len(A) == 4 * (n - 1)
    assert len(B) == (n - 4) ** 2
    assert not any(v in B for u in A for v in neighbors(u, n))

    center = Fraction(n - 1, 2)
    boundary_variance = sum((Fraction(x) - center) ** 2 + (Fraction(y) - center) ** 2
                            for x, y in A) / len(A)
    expected_variance = Fraction(2 * n * n - 4 * n + 3, 6)
    assert boundary_variance == expected_variance

    bridges = {v for v in vertices - A - B
               if any(q in A for q in neighbors(v, n))
               and any(q in B for q in neighbors(v, n))}
    assert len(bridges) == 4 * (n - 4), (n, bridges)

    weights = []
    t_values = []
    for w in bridges:
        S = A | B | {w}
        remainder = vertices - S
        assert len(remainder) == 4 * n - 13
        assert all(any(q in A for q in neighbors(v, n)) for v in remainder)
        # h_G(S)=|remainder|! follows because every remaining point is already attached.
        t_values.append(len(remainder))
        weights.append((len(A), len(B), len(remainder)))
    assert len(weights) == 4 * (n - 4)
    assert len(set(weights)) == 1
    return {"n": n, "boundary_var": str(boundary_variance),
            "equal_weight_states": len(weights), "suffix_size": t_values[0]}


def number_of_births(word, n):
    seen = set()
    births = 0
    for v in word:
        if not any(q in seen for q in neighbors(v, n)):
            births += 1
        seen.add(v)
    return births


def late_merge_ratio(a, b, g, c0):
    n_total = a + b + g
    return Fraction(c0, g) * Fraction(comb(a + b, a), comb(n_total - 1, b - 1))


def audit_n5_concentric_words():
    """Brute-check the two local word families after fixing one A ordering."""
    n = 5
    A_order = [(x, 0) for x in range(n)]
    A_order += [(n - 1, y) for y in range(1, n)]
    A_order += [(x, n - 1) for x in range(n - 2, -1, -1)]
    A_order += [(0, y) for y in range(n - 2, 0, -1)]
    A = set(A_order)
    B = (2, 2)
    bridges = [(1, 2), (2, 1), (3, 2), (2, 3)]
    vertices = {(x, y) for x in range(n) for y in range(n)}
    count_c = 0
    for w in bridges:
        gap = sorted(vertices - A - {B, w})
        for insertion in range(len(A_order) + 1):
            prefix = A_order[:insertion] + [B] + A_order[insertion:] + [w]
            for tail in permutations(gap):
                assert number_of_births(prefix + list(tail), n) == 2
                count_c += 1
    assert count_c == 4 * 17 * factorial(7)

    gap_all = sorted(vertices - A - {B})
    count_d = 0
    for tail in permutations(gap_all):
        assert number_of_births([B] + A_order + list(tail), n) == 2
        count_d += 1
    assert count_d == factorial(8)
    assert late_merge_ratio(16, 1, 8, 4) == Fraction(count_c, count_d)
    return {"n": 5, "C_over_aA": count_c, "D_over_aA": count_d,
            "ratio_C_over_D": str(Fraction(count_c, count_d))}


def audit_rectangle_ratio():
    rows = []
    for n, k in ((12, 5), (20, 9), (32, 15)):
        a, b, g = k * n, (n - k - 1) * n, n
        ratio = late_merge_ratio(a, b, g, g)
        assert 0 < ratio < 1
        rows.append({"n": n, "k": k, "exact_ratio": str(ratio)})
    return rows


if __name__ == "__main__":
    for n in (3, 4, 5, 8, 16):
        if n <= 8:
            brute_path_counts(n)
        print(check(n))
    for n in (5, 6, 8, 16):
        print(concentric_state(n))
    print(audit_n5_concentric_words())
    for row in audit_rectangle_ratio():
        print(row)
