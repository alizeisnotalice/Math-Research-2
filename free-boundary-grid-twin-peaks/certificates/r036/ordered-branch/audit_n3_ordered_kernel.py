#!/usr/bin/env python3
"""Complete P3-square audit of the ordered first-merge kernel."""

from fractions import Fraction
from itertools import permutations
import json
import platform
import sys


N = 3
VERTICES = tuple((x, y) for x in range(N) for y in range(N))


def adjacent(x, y):
    return abs(x[0] - y[0]) + abs(x[1] - y[1]) == 1


NEIGHBOURS = {x: tuple(y for y in VERTICES if adjacent(x, y)) for x in VERTICES}


def first_merge_data(word):
    """Return (a,b,q,rho,e,L) for a K=2 word, or None."""
    owner = {}
    births = []
    roots = []
    merged = False
    result = None
    for rank, v in enumerate(word, start=1):
        touched = {owner[u] for u in NEIGHBOURS[v] if u in owner}
        if not touched:
            births.append(v)
            if len(births) > 2:
                return None
            comp = len(births) - 1
            owner[v] = comp
            roots.append(v)
            continue
        if len(touched) == 1:
            owner[v] = next(iter(touched))
            continue
        if len(touched) != 2 or merged:
            raise AssertionError("unexpected second merge before component relabel")
        if len(births) != 2:
            raise AssertionError("merge before two births")
        sizes = [sum(c == j for c in owner.values()) for j in (0, 1)]
        a, b = sizes[0], sizes[1]  # component 0 contains word[0], the high root
        r1, r2 = roots
        first_low_neighbour_rank = next(
            j for j, z in enumerate(word, start=1)
            if j >= 3 and z in NEIGHBOURS[r2]
        )
        rho = abs(r1[0] - r2[0]) + abs(r1[1] - r2[1])
        e = (r1[0] - r2[0]) ** 2 + (r1[1] - r2[1]) ** 2
        result = (a, b, a + b, rho, e, first_low_neighbour_rank - 2)
        for x in list(owner):
            owner[x] = 0
        owner[v] = 0
        merged = True
    if len(births) != 2 or not merged:
        return None
    return result


def main():
    a2 = 0
    b_total = 0
    full = {}
    top = {}
    lsum = {}
    kernel_total = Fraction()
    reciprocal_high_total = Fraction()
    high_small_rho = {k: 0 for k in range(1, 8)}
    theta_event = {d: 0 for d in (Fraction(1, 4), Fraction(1, 3), Fraction(2, 5))}
    for word in permutations(VERTICES):
        data = first_merge_data(word)
        if data is None:
            continue
        a2 += 1
        a, b, q, rho, e, ell = data
        full[(a, b)] = full.get((a, b), 0) + 1
        reciprocal_high_total += Fraction(q, a)
        for k in high_small_rho:
            if a <= k:
                high_small_rho[k] += rho
        for d in theta_event:
            if Fraction(a, q) <= d:
                theta_event[d] += 1
        # top-two means precisely that the second point is the second birth.
        if word[1] in NEIGHBOURS[word[0]]:
            continue
        b_total += 1
        top[(a, b)] = top.get((a, b), 0) + 1
        lsum[(a, b)] = lsum.get((a, b), 0) + ell
        kernel_total += Fraction(q - 1, b)

    assert a2 == 143808
    assert b_total == 81984
    for (a, b), count in top.items():
        q = a + b
        assert Fraction(lsum[(a, b)], count) == Fraction(q - 1, b)
        assert full[(a, b)] == Fraction(q - 1, b) * count
    assert kernel_total == a2
    assert reciprocal_high_total == 2 * a2
    for d, count in theta_event.items():
        assert Fraction(count, a2) <= d
    for k, rho_sum in high_small_rho.items():
        assert Fraction(rho_sum, a2) <= 2 * k

    print(json.dumps({
        "status": "PASS",
        "python": sys.version.split()[0],
        "implementation": platform.python_implementation(),
        "A2": a2,
        "B": b_total,
        "sum_ordered_kernel_over_B": str(kernel_total),
        "sum_q_over_high_branch_over_A2": str(reciprocal_high_total),
        "top_fibres": {f"{a},{b}": top[(a, b)] for a, b in sorted(top)},
        "theta_event_counts": {str(d): count for d, count in theta_event.items()},
        "high_small_rho_sums": high_small_rho,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
