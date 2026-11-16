#!/usr/bin/env python3
"""Exact audits for the first-merge conditional mean of L on the 3x3 grid."""
from collections import defaultdict
from fractions import Fraction
from itertools import combinations, permutations
import json
from pathlib import Path

N = 9
HERE = Path(__file__).resolve().parent


def xy(v):
    return v % 3, v // 3


ADJ = [set() for _ in range(N)]
for v in range(N):
    x, y = xy(v)
    for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        xx, yy = x + dx, y + dy
        if 0 <= xx < 3 and 0 <= yy < 3:
            ADJ[v].add(yy * 3 + xx)


def births(word):
    seen, out = set(), []
    for rank, v in enumerate(word, 1):
        if not (ADJ[v] & seen):
            out.append((rank, v))
        seen.add(v)
    return out


def components(vertices):
    left, out = set(vertices), []
    while left:
        todo, comp = [left.pop()], set()
        while todo:
            v = todo.pop()
            comp.add(v)
            new = ADJ[v] & left
            left -= new
            todo.extend(new)
        out.append(frozenset(comp))
    return out


def first_merge_data(word):
    """Return the complete fixed first-merge fibre datum and its direction."""
    b = births(word)
    assert len(b) == 2
    tau = b[1][0]
    for r in range(tau + 1, N + 1):
        v, old = word[r - 1], set(word[:r - 1])
        cc = components(old)
        touched = sum(bool(ADJ[v] & c) for c in cc)
        if touched >= 2:
            w, S = v, frozenset(word[:r])
            cs = components(set(S) - {w})
            assert len(cs) == 2
            A, B = sorted(cs, key=lambda c: min(c))
            seqA = tuple(z for z in word[:r - 1] if z in A)
            seqB = tuple(z for z in word[:r - 1] if z in B)
            suffix = tuple(word[r:])
            key = (tuple(sorted(S)), w, seqA, seqB, suffix)
            direction = "A" if word[0] in A else "B"
            return key, len(A), len(B), direction
    raise AssertionError("two components never merged")


def e2(x, y):
    xx, xy0 = xy(x)
    yx, yy = xy(y)
    return (xx - yx) ** 2 + (xy0 - yy) ** 2


def L_tau2(word):
    assert births(word)[1][0] == 2
    y = word[1]
    j = next(i for i, v in enumerate(word, 1) if i >= 3 and v in ADJ[y])
    return j - 2


def binary_mean(first, second, a, b):
    """Mean L after fixed first/second roots, by enumerating all interleavings."""
    # The residual A/B items retain their internal order; their choices of A
    # slots parameterize every such interleaving exactly once.
    q, total = a + b, Fraction()
    cnt = 0
    for Aslots in combinations(range(q - 2), a - 1):
        letters = ["B"] * (q - 2)
        for i in Aslots:
            letters[i] = "A"
        target = second
        if target not in letters:
            # b=1 in AB (or a=1 in BA): the merge vertex is the first
            # neighbour of the second root, at rank q+1.
            L = q - 1
        else:
            L = letters.index(target) + 1
        total += L
        cnt += 1
    return total / cnt


def main():
    # Pure interleaving audit, all a,b <= 8, with exact rational means.
    checked_pairs = 0
    for a in range(1, 9):
        for b in range(1, 9):
            q = a + b
            mean_ab = binary_mean("A", "B", a, b)
            mean_ba = binary_mean("B", "A", a, b)
            assert mean_ab == Fraction(q - 1, b)
            assert mean_ba == Fraction(q - 1, a)
            checked_pairs += 1

    # Full words are grouped by the complete first-merge fibre datum.  On the
    # tau=2 subcollection, L must recover every full K=2 word, direction by
    # direction, while 1/theta is the root-symmetric Rao--Blackwell weight.
    full = defaultdict(lambda: [0, 0])
    tau2_L = defaultdict(lambda: [0, 0])
    tau2_theta = defaultdict(lambda: [Fraction(), Fraction()])
    total_D = 0
    for word in permutations(range(N)):
        b = births(word)
        if len(b) != 2:
            continue
        key, a, bb, direction = first_merge_data(word)
        d = e2(b[0][1], b[1][1])
        full[(key, direction)][0] += 1
        full[(key, direction)][1] += d
        total_D += d
        if b[1][0] != 2:
            continue
        L = L_tau2(word)
        tau2_L[(key, direction)][0] += L
        tau2_L[(key, direction)][1] += L * d
        theta_inv = Fraction((a + bb) * (a + bb - 1), 2 * a * bb)
        tau2_theta[key][0] += theta_inv
        tau2_theta[key][1] += theta_inv * d

    assert set(full) == set(tau2_L)
    assert all(full[k] == tau2_L[k] for k in full)
    # Aggregating directions is exactly where the symmetric 1/theta weight
    # applies.  Check it fibre by fibre, not only after the global sum.
    by_key_full = defaultdict(lambda: [0, 0])
    for (key, _direction), val in full.items():
        by_key_full[key][0] += val[0]
        by_key_full[key][1] += val[1]
    assert set(by_key_full) == set(tau2_theta)
    assert all(tuple(by_key_full[k]) == tuple(tau2_theta[k]) for k in by_key_full)

    output = {
        "n": 3,
        "method": "all 9! words; complete first-merge fibre keys and exact Fraction interleavings",
        "small_interleaving_pairs_checked": checked_pairs,
        "complete_fibres": len(by_key_full),
        "direct_A2": sum(v[0] for v in full.values()),
        "directionwise_L_weighted_A2": sum(v[0] for v in tau2_L.values()),
        "direct_euclidean_square_total": total_D,
        "directionwise_L_weighted_euclidean_square_total": sum(v[1] for v in tau2_L.values()),
        "symmetric_theta_inverse_A2": str(sum(v[0] for v in tau2_theta.values())),
        "symmetric_theta_inverse_euclidean_square_total": str(sum(v[1] for v in tau2_theta.values())),
        "checks": {
            "all_small_fraction_means": True,
            "every_direction_fibre_L_recovers_full_count": True,
            "every_direction_fibre_L_recovers_full_D": True,
            "every_fibre_theta_inverse_recovers_full_count": True,
            "every_fibre_theta_inverse_recovers_full_D": True,
        },
    }
    assert output["direct_A2"] == 143808
    assert output["direct_euclidean_square_total"] == 633600
    (HERE / "first-merge-rb-n3.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
