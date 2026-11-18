"""Exact small-dimensional checks for the crowded-plane criterion."""
from collections import Counter, defaultdict
from itertools import combinations
import json
import platform


def walsh(values):
    out = list(values)
    width = 1
    while width < len(out):
        for base in range(0, len(out), 2 * width):
            for i in range(base, base + width):
                a, b = out[i], out[i + width]
                out[i], out[i + width] = a + b, a - b
        width *= 2
    return out


def geometry(points):
    counts = Counter(x ^ y for x, y in combinations(points, 2))
    planes = []
    for x, y, z in combinations(sorted(points), 3):
        w = x ^ y ^ z
        if w > z and w in points:
            mask = (1 << x) | (1 << y) | (1 << z) | (1 << w)
            crowded = min(counts[x ^ y], counts[x ^ z], counts[x ^ w]) >= 3
            planes.append((mask, crowded))
    return counts, planes


def rank(points):
    pivots = {}
    for x in points:
        while x:
            i = x.bit_length()
            if i not in pivots:
                pivots[i] = x
                break
            x ^= pivots[i]
    return len(pivots)


def audit(n):
    size = 1 << n
    # Modding out by a global sign leaves the sign parity of four points unchanged.
    saved = {}
    types = Counter()
    least_positive = defaultdict(lambda: 10**9)
    for mask in range(1 << (size - 1)):
        c = walsh(1 - 2 * (mask >> x & 1) for x in range(size))
        T = frozenset(x for x, z in enumerate(c) if z)
        if T not in saved:
            saved[T] = geometry(T)
        counts, planes = saved[T]
        neg = sum(1 << x for x in T if c[x] < 0)
        positive = [m for m, crowded in planes if (m & neg).bit_count() % 2 == 0]
        assert all(crowded or (m & neg).bit_count() % 2 == 1 for m, crowded in planes)
        if len(T) > 4:
            assert all(any(counts[x ^ z] >= 3 for z in T if x != z) for x in T)
            assert len(T) <= 10 * len(positive)
            assert rank(T) <= 6 * len(positive)
        else:
            assert len(T) in (1, 4)
            if len(T) == 4:
                assert len(planes) == 1 and all(abs(c[x]) == size // 2 for x in T)
        if not any(crowded for _, crowded in planes):
            assert len(T) in (1, 4)
        types[len(T)] += 2
        least_positive[len(T)] = min(least_positive[len(T)], len(positive))
    return {"boolean_tables": sum(types.values()), "distinct_supports": len(saved),
            "by_support_size": dict(sorted(types.items())),
            "min_positive_planes_per_size": dict(sorted(least_positive.items()))}


def independent_blocks(blocks):
    S = {0}
    for i in range(blocks):
        unit = 1 << (3 * i)
        S.update(unit * c for c in (1, 3, 5, 7))
    # Avoid a bitmask indexed by high-dimensional frequencies.
    by_delta = defaultdict(list)
    for a, b in combinations(S, 2):
        by_delta[a ^ b].append((a, b))
    assert max(map(len, by_delta.values())) == 2
    return {"m": blocks, "q": 4 * blocks**2, "r": rank(S),
            "max_unordered_pair_multiplicity": 2}


def crowded_examples(groups, width):
    # Separate packets of independent parallel pairs with unrelated crowded directions.
    S = {0}
    for g in range(groups):
        shift = g * (width + 1)
        lam = 1 << shift
        for i in range(width):
            t = 1 << (shift + i + 1)
            S.update((t, t ^ lam))
    counts = Counter(x ^ y for x, y in combinations(S, 2))
    D = {d for d, v in counts.items() if v >= 3}
    assert len(D) == groups and all(counts[d] == width for d in D)
    assert all((d ^ e) not in D for d, e in combinations(D, 2))
    return {"groups": groups, "pairs_per_crowded_direction": width, "rank_of_crowded_directions": rank(D)}


print(json.dumps({"python": platform.python_version(),
                  "complete_census": {str(n): audit(n) for n in range(1, 5)},
                  "equality_family": [independent_blocks(1 << j) for j in range(5)],
                  "high_multiplicity_sum_free": crowded_examples(4, 10)}, indent=2))
