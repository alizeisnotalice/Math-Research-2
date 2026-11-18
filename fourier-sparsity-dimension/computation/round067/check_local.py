"""All candidate label sets in dimension four, using exact bit masks."""
from collections import Counter
from functools import cache
from itertools import combinations
import json
from check_peeling import support

SIZE = 16
PAIRS = [[(1 << x) | (1 << (x ^ d)) for x in range(SIZE) if x < (x ^ d)]
         for d in range(SIZE)]


@cache
def counts(mask):
    return tuple(sum(mask & p == p for p in terms) for terms in PAIRS)


@cache
def admits_crowded_plane(mask):
    heavy = [d for d, count in enumerate(counts(mask)) if count >= 3]
    return any(mask & plane == plane
               for a, b in combinations(heavy, 2) if a ^ b in heavy
               for x in range(SIZE)
               for plane in [(1 << x) | (1 << (x ^ a)) |
                             (1 << (x ^ b)) | (1 << (x ^ a ^ b))])


def local_core(mask):
    steps = 0
    while mask:
        retained = 0
        m = counts(mask)
        for a in range(SIZE):
            if not mask >> a & 1:
                continue
            neighborhood = sum(1 << b for b in range(SIZE) if mask >> b & 1
                               and (a == b or m[a ^ b] >= 2))
            local = counts(neighborhood)
            if any(local[a ^ b] >= 3 for b in range(SIZE)
                   if neighborhood >> b & 1) and admits_crowded_plane(neighborhood):
                retained |= 1 << a
        if retained == mask:
            break
        mask = retained
        steps += 1
    return mask, steps


if __name__ == '__main__':
    unions = [0] * (1 << SIZE)
    for table in range(1 << SIZE):
        labels = support(table, SIZE)
        if len(labels) > 4:
            mask = sum(1 << x for x in labels)
            unions[mask] = mask
    steps_hist = Counter()
    strict = example = 0
    for mask in range(1 << SIZE):
        for x in range(SIZE):
            if mask >> x & 1:
                unions[mask] |= unions[mask ^ (1 << x)]
        kept, steps = local_core(mask)
        assert not unions[mask] & ~kept
        basic = sum(1 << a for a in range(SIZE) if mask >> a & 1 and
                    any(counts(mask)[a ^ b] >= 3 for b in range(SIZE) if mask >> b & 1))
        assert not kept & ~basic
        if basic != kept:
            strict += 1
            if not example:
                example = dict(mask=mask, basic=basic, localized=kept)
        steps_hist[steps] += 1
    print(json.dumps(dict(ambient=4, tested=len(unions), strict_improvements=strict,
                          example=example, rounds=dict(sorted(steps_hist.items()))), indent=2))
