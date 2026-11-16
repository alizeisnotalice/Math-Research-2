#!/usr/bin/env python3
"""Independent exhaustive n=3 one-peak oracle for ordered geodesic prefixes."""
import itertools
import json

N_SIDE = 3
N = N_SIDE * N_SIDE
MAX_D = 2 * N_SIDE - 1

def manhattan(a, b):
    return abs(a % N_SIDE - b % N_SIDE) + abs(a // N_SIDE - b // N_SIDE)

def neighbours(v):
    x, y = v % N_SIDE, v // N_SIDE
    for dx, dy in ((-1, 0), (1, 0), (0, -1), (0, 1)):
        xx, yy = x + dx, y + dy
        if 0 <= xx < N_SIDE and 0 <= yy < N_SIDE:
            yield yy * N_SIDE + xx

def geodesic_prefix(word, d):
    return (all(manhattan(word[k - 1], word[k]) == 1 for k in range(1, d))
            and manhattan(word[0], word[d - 1]) == d - 1)

counts = [0] * MAX_D
one_peak_words = 0
for word in itertools.permutations(range(N)):
    position = [0] * N
    for rank, vertex in enumerate(word):
        position[vertex] = rank
    births = sum(not any(position[w] < position[v] for w in neighbours(v))
                 for v in range(N))
    if births != 1:
        continue
    one_peak_words += 1
    for d in range(1, MAX_D + 1):
        counts[d - 1] += geodesic_prefix(word, d)

assert one_peak_words == 29568, one_peak_words
assert counts == [29568, 29568, 14784, 4000, 448], counts
print(json.dumps({
    "kind": "independent-n3-onepeak-geodesic-oracle",
    "words": one_peak_words,
    "d": list(range(1, MAX_D + 1)),
    "counts": counts,
    "probabilities": [x / one_peak_words for x in counts],
}, indent=2))
