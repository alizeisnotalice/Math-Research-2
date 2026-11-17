#!/usr/bin/env python3
"""Exhaustive n=3 checks for the two deterministic word maps in R029."""
from __future__ import annotations

import itertools
import json
from collections import Counter
from math import factorial


def graph(n: int):
    adj = [set() for _ in range(n * n)]
    for y in range(n):
        for x in range(n):
            v = y * n + x
            for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                xx, yy = x + dx, y + dy
                if 0 <= xx < n and 0 <= yy < n:
                    adj[v].add(yy * n + xx)
    return adj


def peaks(word, adj):
    pos = {x: i for i, x in enumerate(word)}
    return [x for x in word if all(pos[x] < pos[y] for y in adj[x] if y in pos)]


def corner_delete_map(word):
    # p=(0,0), D={0,1,3}; T is the connected path 2-5-4-7-6.
    D, T = {0, 1, 3}, {2, 5, 4, 7, 6}
    C, root = D | T, word[0]
    outside = [x for x in word if x not in C]
    if root in C:
        return tuple([2, 5, 4, 7, 6] + outside)
    c = next(x for x in word if x in C)
    assert c in T
    path, i = [2, 5, 4, 7, 6], [2, 5, 4, 7, 6].index(c)
    block = [c]
    for d in range(1, len(path)):
        if i - d >= 0:
            block.append(path[i - d])
        if i + d < len(path):
            block.append(path[i + d])
    j = word.index(c)
    return tuple([x for x in word[:j] if x not in C] + block +
                 [x for x in word[j + 1:] if x not in C])


def canonical_path(r, q, n):
    x, y, xx, yy = r % n, r // n, q % n, q // n
    ans = [r]
    while x != xx:
        x += 1 if xx > x else -1
        ans.append(y * n + x)
    while y != yy:
        y += 1 if yy > y else -1
        ans.append(y * n + x)
    return ans


def path_promote(word, n):
    path = canonical_path(word[0], word[1], n)
    return tuple(path + [x for x in word if x not in set(path)])


def main():
    n, N = 3, 9
    adj = graph(n)
    H_adj = [a - {0, 1, 3} for a in adj]
    local_images, path_images = Counter(), Counter()
    a1 = b = 0
    for word in itertools.permutations(range(N)):
        pk = peaks(word, adj)
        if len(pk) == 1:
            a1 += 1
            image = corner_delete_map(word)
            assert set(image) == {2, 4, 5, 6, 7, 8}
            assert len(peaks(image, H_adj)) == 1, (word, image)
            local_images[image] += 1
        if len(pk) == 2 and word[1] not in adj[word[0]]:
            b += 1
            image = path_promote(word, n)
            assert len(peaks(image, adj)) == 1, (word, image)
            path_images[image] += 1
    print(json.dumps({
        "status": "passed", "n": n, "A1_words_checked": a1,
        "corner_delete_max_fibre": max(local_images.values()),
        "corner_delete_bound_falling_N_8": factorial(N),
        "B_top2_words_checked": b,
        "path_promote_max_fibre": max(path_images.values()),
        "path_promote_bound_N_power_Dmax": N ** (2 * n - 2),
    }, indent=2))


if __name__ == "__main__":
    main()
