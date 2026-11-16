#!/usr/bin/env python3
"""Exact n=3 tests for two tempting local root-transport maps."""
import json
from collections import Counter
from pathlib import Path


n = 3
N = n * n
corner, neighbour = 0, 1
adj = [set() for _ in range(N)]
for v in range(N):
    x, y = v % n, v // n
    for xx, yy in ((x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1)):
        if 0 <= xx < n and 0 <= yy < n:
            adj[v].add(yy * n + xx)


def connected_word(word):
    seen = set()
    for v in word:
        if seen and not (adj[v] & seen):
            return False
        seen.add(v)
    return True


def words(root):
    out = []
    def walk(prefix, seen):
        if len(prefix) == N:
            out.append(tuple(prefix))
            return
        for v in range(N):
            if v not in seen and adj[v] & seen:
                walk(prefix + [v], seen | {v})
    walk([root], {root})
    return out


corner_words = words(corner)
neighbour_words = words(neighbour)
assert len(corner_words) == 1344 and len(neighbour_words) == 3696

swap_good = 0
stable_images = Counter()
swap_bad_witness = None
for word in corner_words:
    j = word.index(neighbour)
    swapped = list(word)
    swapped[0], swapped[j] = swapped[j], swapped[0]
    if connected_word(swapped) and swapped[0] == neighbour:
        swap_good += 1
    elif swap_bad_witness is None:
        swap_bad_witness = {"input": list(word), "swapped": swapped}
    stable = (neighbour,) + tuple(v for v in word if v != neighbour)
    assert connected_word(stable) and stable[:2] == (neighbour, corner)
    stable_images[stable] += 1

payload = {
    "n": n,
    "corner_words": len(corner_words),
    "neighbour_words": len(neighbour_words),
    "exact_ratio": "11/4",
    "label_position_transposition_good": swap_good,
    "label_position_transposition_bad": len(corner_words) - swap_good,
    "label_position_transposition_bad_witness": swap_bad_witness,
    "stable_move_image_size": len(stable_images),
    "stable_move_max_fibre": max(stable_images.values()),
    "stable_move_fibre_histogram": {str(k): v for k, v in sorted(Counter(stable_images.values()).items())},
}
Path(__file__).with_name("coordinate-unimodality-n3.json").write_text(json.dumps(payload, indent=2) + "\n")
print(json.dumps(payload, indent=2))
