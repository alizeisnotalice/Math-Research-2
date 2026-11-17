#!/usr/bin/env python3
"""Independent finite audit for the corner-deletion connected-word map.

Coordinates are (x,y), 0 <= x,y < n.  A word is read from largest to
smallest label; it represents a one-peak labeling precisely when every entry
after its first has a neighbour earlier in the word.
"""

from collections import Counter, deque
from itertools import permutations
from math import prod


def vertices(n):
    return tuple((x, y) for y in range(n) for x in range(n))


def adjacent(u, v):
    return abs(u[0] - v[0]) + abs(u[1] - v[1]) == 1


def neighbourhood(vs):
    return {v: {w for w in vs if adjacent(v, w)} for v in vs}


def connected_word(word, nbr):
    seen = {word[0]}
    for v in word[1:]:
        if not (nbr[v] & seen):
            return False
        seen.add(v)
    return True


def peaks_in_descending_word(word, nbr):
    """Vertices with no earlier neighbour, hence peaks of the associated labels."""
    seen = set()
    ans = []
    for v in word:
        if not (nbr[v] & seen):
            ans.append(v)
        seen.add(v)
    return tuple(ans)


def rooted_connected_order_count(vs, nbr):
    """Exact subset DP for all connected orders (one-peak labelings)."""
    index = {v: i for i, v in enumerate(vs)}
    full = (1 << len(vs)) - 1
    masks_nbr = [sum(1 << index[w] for w in nbr[v]) for v in vs]
    dp = [0] * (full + 1)
    for i in range(len(vs)):
        dp[1 << i] = 1
    for mask in range(1, full + 1):
        value = dp[mask]
        if not value:
            continue
        for i in range(len(vs)):
            bit = 1 << i
            if not (mask & bit) and (masks_nbr[i] & mask):
                dp[mask | bit] += value
    return dp[full]


def bfs_order_from(start, tset, nbr):
    """A deterministic connected order of T beginning at start."""
    seen = {start}
    queue = deque([start])
    out = []
    while queue:
        v = queue.popleft()
        out.append(v)
        for w in sorted(nbr[v] & tset):
            if w not in seen:
                seen.add(w)
                queue.append(w)
    assert set(out) == set(tset)
    return tuple(out)


def map_word(word, cset, tset, uset, nbr):
    """The proposed deletion/reinsertion map G-words -> H-words."""
    tau = bfs_order_from((2, 0), tset, nbr)
    if word[0] in cset:
        return tau + tuple(v for v in word if v in uset)
    first_c_index = next(i for i, v in enumerate(word) if v in cset)
    c = word[first_c_index]
    assert c in tset, (word, c)
    tau_c = bfs_order_from(c, tset, nbr)
    # Before c all entries are U.  Delete all original C vertices and put T
    # at the former location of c, while retaining the U relative order.
    prefix_u = tuple(v for v in word[:first_c_index] if v in uset)
    suffix_u = tuple(v for v in word[first_c_index:] if v in uset)
    return prefix_u + tau_c + suffix_u


def exhaustive_n3():
    n = 3
    vs = vertices(n)
    nbr_g = neighbourhood(vs)
    dset = {(0, 0), (1, 0), (0, 1)}
    tset = {(2, 0), (2, 1), (1, 1), (1, 2), (0, 2)}
    cset = dset | tset
    uset = set(vs) - cset
    hvs = tuple(v for v in vs if v not in dset)
    nbr_h = neighbourhood(hvs)
    one_peak = 0
    all_permutations = 0
    images = Counter()
    root_classes = Counter()
    for word in permutations(vs):
        all_permutations += 1
        peaks = peaks_in_descending_word(word, nbr_g)
        if len(peaks) != 1:
            continue
        # Audit the equivalence used by the map, not merely a recurrence.
        assert connected_word(word, nbr_g)
        one_peak += 1
        root_classes["C" if word[0] in cset else "U"] += 1
        out = map_word(word, cset, tset, uset, nbr_g)
        assert set(out) == set(hvs) and len(out) == len(hvs)
        assert connected_word(out, nbr_h), (word, out)
        assert len(peaks_in_descending_word(out, nbr_h)) == 1, (word, out)
        images[out] += 1
    h_enum = sum(1 for word in permutations(hvs) if len(peaks_in_descending_word(word, nbr_h)) == 1)
    assert all_permutations == 362880
    assert one_peak == rooted_connected_order_count(vs, nbr_g)
    assert h_enum == rooted_connected_order_count(hvs, nbr_h)
    return {
        "all_label_orders_examined": all_permutations,
        "A1_G": one_peak,
        "A1_H": h_enum,
        "root_in_C": root_classes["C"],
        "root_in_U": root_classes["U"],
        "distinct_H_images": len(images),
        "max_observed_fibre": max(images.values()),
        "falling_9_8": prod(range(2, 10)),
        "inequality": one_peak <= prod(range(2, 10)) * h_enum,
    }


def subset_dp_n4():
    n = 4
    vs = vertices(n)
    dset = {(0, 0), (1, 0), (0, 1)}
    hvs = tuple(v for v in vs if v not in dset)
    a_g = rooted_connected_order_count(vs, neighbourhood(vs))
    a_h = rooted_connected_order_count(hvs, neighbourhood(hvs))
    factor = prod(range(n * n - 7, n * n + 1))
    return {
        "A1_G": a_g,
        "A1_H": a_h,
        "falling_16_8": factor,
        "rhs": factor * a_h,
        "inequality": a_g <= factor * a_h,
    }


def main():
    print("corner-deletion map exact audit")
    print("n=3 exhaustive permutations:")
    for key, value in exhaustive_n3().items():
        print(f"  {key}: {value}")
    print("n=4 subset dynamic programming:")
    for key, value in subset_dp_n4().items():
        print(f"  {key}: {value}")


if __name__ == "__main__":
    main()
