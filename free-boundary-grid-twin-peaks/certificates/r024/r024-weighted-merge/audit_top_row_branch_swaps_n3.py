#!/usr/bin/env python3
"""Full 3x3 audit of the four reversible top-row branch moves."""
from itertools import permutations
import json
from pathlib import Path

n = 3
N = n * n
TOP = {0, 1, 2}
ROW1 = {3, 4, 5}


def neighbors(v):
    x, y = v % n, v // n
    return {Y * n + X for X, Y in ((x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1))
            if 0 <= X < n and 0 <= Y < n}


ADJ = [neighbors(v) for v in range(N)]


def peaks(word):
    rank = [0] * N
    for i, v in enumerate(word):
        rank[v] = i
    return tuple(v for v in range(N) if all(rank[v] < rank[u] for u in ADJ[v]))


def first_merge(word):
    seen = set()
    components = []
    for r, v in enumerate(word, 1):
        hit = [i for i, c in enumerate(components) if ADJ[v] & c]
        if not hit:
            components.append({v})
        else:
            merged = {v}
            for i in reversed(hit):
                merged |= components.pop(i)
            components.append(merged)
        if len(components) == 1 and r > 1:
            # The first time two pre-existing components coalesce occurs only
            # when the preceding state had exactly two components.
            prev = set(word[:r - 1])
            # Count its components directly.
            todo = set(prev)
            ccount = 0
            while todo:
                ccount += 1
                stack = {next(iter(todo))}
                while stack:
                    u = stack.pop()
                    if u not in todo:
                        continue
                    todo.remove(u)
                    stack |= ADJ[u] & todo
            if ccount == 2:
                return r, frozenset(word[:r]), v
    return None


def swap_positions(word, cols):
    out = list(word)
    pos = [0] * N
    for i, v in enumerate(out):
        pos[v] = i
    for x in cols:
        u, v = x, n + x
        iu, iv = pos[u], pos[v]
        out[iu], out[iv] = out[iv], out[iu]
    return tuple(out)


def swap_rows_on_columns(word, row, cols):
    out = list(word)
    pos = [0] * N
    for i, v in enumerate(out):
        pos[v] = i
    for x in cols:
        u, v = row * n + x, (row + 1) * n + x
        iu, iv = pos[u], pos[v]
        out[iu], out[iv] = out[iv], out[iu]
    return tuple(out)


source = []
for word in permutations(range(N)):
    if len(peaks(word)) != 2:
        continue
    event = first_merge(word)
    if event and event[0] == n and event[1] == TOP:
        source.append(word)

assert len(source) == 312  # Exact F_3 computed independently by h recursion.
images = [set() for _ in range(4)]
for word in source:
    event = first_merge(word)
    j = event[2]  # bridge column in the top row
    masks = [set(), set(range(j)), set(range(j + 1, n)), set(range(n))]
    for q, cols in enumerate(masks):
        image = swap_positions(word, cols)
        assert len(peaks(image)) == 2
        images[q].add(image)
assert all(len(x) == len(source) for x in images)
assert all(images[i].isdisjoint(images[j]) for i in range(4) for j in range(i))

def peaks_grid(word, side):
    rank = [0] * (side * side)
    for i, v in enumerate(word):
        rank[v] = i
    ans = []
    for v in range(side * side):
        x, y = v % side, v // side
        ns = [Y * side + X for X, Y in ((x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1))
              if 0 <= X < side and 0 <= Y < side]
        if all(rank[v] < rank[u] for u in ns):
            ans.append(v)
    return ans

# The n=3 oracle happens to tolerate one extra step.  The following explicit
# n=4 valid top-row arm word shows that a second partial move is not a
# uniform construction in n: it creates a third peak.
source_n4 = tuple([0, 1, 3, 2] + list(range(4, 16)))
once_n4 = list(source_n4)
for u, v in ((0, 4), (1, 5)):
    iu, iv = once_n4.index(u), once_n4.index(v)
    once_n4[iu], once_n4[iv] = once_n4[iv], once_n4[iu]
twice_n4 = list(once_n4)
for u, v in ((4, 8), (5, 9)):
    iu, iv = twice_n4.index(u), twice_n4.index(v)
    twice_n4[iu], twice_n4[iv] = twice_n4[iv], twice_n4[iu]
assert len(peaks_grid(source_n4, 4)) == 2
assert len(peaks_grid(once_n4, 4)) == 2
assert peaks_grid(twice_n4, 4) == [0, 3, 8]
repeat_failure = {"n": 4, "source": list(source_n4), "after_first_move": once_n4,
                  "after_second_move": twice_n4, "peaks_after_second_move": [0, 3, 8]}

out = {"n": n, "source_F_n": len(source), "image_sizes": [len(x) for x in images],
       "pairwise_disjoint": True, "all_images_have_exactly_two_peaks": True,
       "bound_verified_on_oracle": "4F_3 <= A2(3)",
       "second_partial_move_failure": repeat_failure}
Path(__file__).with_name("top_row_branch_swaps_n3.json").write_text(json.dumps(out, indent=2) + "\n")
print(json.dumps(out, indent=2))
