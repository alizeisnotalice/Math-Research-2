#!/usr/bin/env python3
"""Exact subset-DP audit of canonical detachable-branch images on 4x4."""
from collections import defaultdict
from math import comb
import json
from pathlib import Path
import sys

n = int(sys.argv[1]) if len(sys.argv) == 2 else 4
N = n * n
FULL = (1 << N) - 1


def make_adj():
    out = [0] * N
    for y in range(n):
        for x in range(n):
            v = y * n + x
            for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                xx, yy = x + dx, y + dy
                if 0 <= xx < n and 0 <= yy < n:
                    out[v] |= 1 << (yy * n + xx)
    return out


ADJ = make_adj()


def bits(mask):
    while mask:
        low = mask & -mask
        yield low.bit_length() - 1
        mask ^= low


def comps(mask):
    result = []
    unseen = mask
    while unseen:
        start = unseen & -unseen
        unseen ^= start
        component, todo = 0, start
        while todo:
            low = todo & -todo
            todo ^= low
            vertex = low.bit_length() - 1
            component |= low
            fresh = ADJ[vertex] & unseen
            todo |= fresh
            unseen ^= fresh
        result.append(component)
    return result


def main():
    roots = [[0] * N for _ in range(FULL + 1)]
    for size in range(1, N + 1):
        for mask in range(1, FULL + 1):
            if mask.bit_count() != size:
                continue
            if size == 1:
                roots[mask][(mask & -mask).bit_length() - 1] = 1
                continue
            value = roots[mask]
            for v in bits(mask):
                previous = mask ^ (1 << v)
                if ADJ[v] & previous:
                    for u, count in enumerate(roots[previous]):
                        value[u] += count

    h = [0] * (FULL + 1)
    h[FULL] = 1
    for size in range(N - 1, 0, -1):
        for mask in range(1, FULL + 1):
            if mask.bit_count() == size:
                h[mask] = sum(h[mask | (1 << v)] for v in bits(FULL ^ mask) if ADJ[v] & mask)

    a1, a2, beta_positive, image = sum(roots[FULL]), 0, 0, 0
    rank_total, rank_image = defaultdict(int), defaultdict(int)
    size_total, size_image = defaultdict(int), defaultdict(int)
    for mask in range(1, FULL + 1):
        if mask.bit_count() < 3 or not h[mask]:
            continue
        for w in bits(mask):
            pair = comps(mask ^ (1 << w))
            if len(pair) != 2 or not all(ADJ[w] & component for component in pair):
                continue
            A, B = sorted(pair, key=lambda component: (component & -component).bit_length())
            a, b = A.bit_count(), B.bit_count()
            aA, aB = sum(roots[A]), sum(roots[B])
            common = h[mask] * comb(a + b, a)
            weight = common * aA * aB
            a2 += weight
            rank_total[mask.bit_count()] += weight
            size_total[f"{a},{b}"] += weight
            if a == b == 1:
                beta_positive += weight
            elif a == 1:
                beta_positive += h[mask] * aB * b
            elif b == 1:
                beta_positive += h[mask] * aA * a
            # This is exactly the root-B~w condition, not merely B∩N(w) nonempty.
            root_B_near_w = sum(roots[B][v] for v in bits(B & ADJ[w]))
            selected = common * aA * root_B_near_w
            image += selected
            rank_image[mask.bit_count()] += selected
            size_image[f"{a},{b}"] += selected

    expected_values = {
        3: (29568, 143808, 79360),
        4: (35704394880, 779299174080, 266753826624),
    }
    expected = expected_values[n]
    if (a1, a2, beta_positive) != expected:
        raise AssertionError({"computed": (a1, a2, beta_positive), "expected": expected})
    data = {
        "grid": "4x4", "A1": a1, "A2": a2,
        "beta_positive_count": beta_positive, "beta_zero_count": a2 - beta_positive,
        "canonical_detachable_image_count": image,
        "image_over_A1": {"numerator": image, "denominator": a1},
        "image_over_A2": {"numerator": image, "denominator": a2},
        "by_first_merge_rank": {str(k): {"total_W": rank_total[k], "image": rank_image[k]} for k in sorted(rank_total)},
        "by_canonical_branch_sizes": {k: {"total_W": size_total[k], "image": size_image[k]} for k in sorted(size_total, key=lambda x: tuple(map(int, x.split(','))))},
        "status": "passed",
    }
    Path(f"first-merge-root-n{n}.json").write_text(json.dumps(data, indent=2) + "\n")
    print(json.dumps({k: data[k] for k in ("A1", "A2", "beta_positive_count", "beta_zero_count", "canonical_detachable_image_count", "image_over_A1", "image_over_A2", "status")}, indent=2))


if __name__ == "__main__":
    main()
