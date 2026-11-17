#!/usr/bin/env python3
"""Exact local-move audit on P3 square; no asymptotic claim."""

import json
import platform
import sys


N = 3
VERTICES = [(x, y) for x in range(N) for y in range(N)]


def adjacent(a, b):
    return abs(a[0] - b[0]) + abs(a[1] - b[1]) == 1


def components(seen):
    unseen = set(seen)
    answer = []
    while unseen:
        seed = unseen.pop()
        comp = {seed}
        stack = [seed]
        while stack:
            x = stack.pop()
            for y in list(unseen):
                if adjacent(x, y):
                    unseen.remove(y)
                    comp.add(y)
                    stack.append(y)
        answer.append(comp)
    return answer


def trace(word):
    seen = []
    births = []
    first_merge = None
    low_singleton_at_merge = None
    for rank, v in enumerate(word, start=1):
        old_components = components(seen) if seen else []
        touched = [c for c in old_components if any(adjacent(v, u) for u in c)]
        if not touched:
            births.append((rank, v))
        if len(touched) >= 2 and first_merge is None:
            first_merge = rank
            high_root, low_root = word[0], word[1]
            high_comp = next(c for c in old_components if high_root in c)
            low_comp = next(c for c in old_components if low_root in c)
            low_singleton_at_merge = len(low_comp) == 1
        seen.append(v)
    return {
        "births": births,
        "first_merge_rank": first_merge,
        "low_singleton_at_first_merge": low_singleton_at_merge,
    }


def main():
    # v is a corner low root; w is its earliest neighbour and first merger.
    # y is the only other neighbour of v.  Moving y just before w is the
    # natural one-step attempt to grow the low branch, but y also touches A.
    u, v = (2, 2), (0, 0)
    w, y = (1, 0), (0, 1)
    original = [u, v, (2, 1), (1, 2), (2, 0), (1, 1), (0, 2), w, y]
    swapped = [u, v, (2, 1), (1, 2), (2, 0), (1, 1), (0, 2), y, w]
    assert set(original) == set(VERTICES) and len(original) == len(VERTICES)
    assert set(swapped) == set(VERTICES) and len(swapped) == len(VERTICES)
    assert set(x for x in VERTICES if adjacent(v, x)) == {w, y}
    result = {"original": trace(original), "swap_y_before_w": trace(swapped)}
    assert len(result["original"]["births"]) == 2
    assert len(result["swap_y_before_w"]["births"]) == 2
    assert result["original"]["first_merge_rank"] == 8
    assert result["swap_y_before_w"]["first_merge_rank"] == 8
    assert result["original"]["low_singleton_at_first_merge"]
    assert result["swap_y_before_w"]["low_singleton_at_first_merge"]
    print(json.dumps({
        "status": "PASS",
        "python": sys.version.split()[0],
        "implementation": platform.python_implementation(),
        "word": original,
        "alternative_neighbour": y,
        "result": result,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
