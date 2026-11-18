from collections import Counter, deque
from fractions import Fraction as F
from itertools import permutations, product
from random import Random


def grid(n):
    points = [(x, y) for y in range(n) for x in range(n)]
    ids = {p: i for i, p in enumerate(points)}
    neighbors = [tuple(ids[q] for q in ((x-1,y),(x+1,y),(x,y-1),(x,y+1)) if q in ids) for x, y in points]
    masks = [sum(1 << j for j in a) for a in neighbors]
    return points, ids, neighbors, masks


def distance(a, b, points):
    return sum(abs(s-t) for s, t in zip(points[a], points[b]))


def peaks(word, masks):
    seen = 0
    out = []
    for v in word:
        if not (masks[v] & seen):
            out.append(v)
        seen |= 1 << v
    return out


def first_merge(word, masks):
    components = []
    for v in word:
        touched = [i for i, mask in enumerate(components) if masks[v] & mask]
        if len(touched) == 2:
            return tuple(components[i] for i in touched)
        if touched:
            components[touched[0]] |= 1 << v
        else:
            components.append(1 << v)
    raise AssertionError("two-peak word never merged")


def check_word(word, points, neighbors, masks):
    roots = peaks(word, masks)
    if len(roots) != 2:
        return None
    p, q = roots
    d = distance(p, q, points)
    if d < 3:
        return None
    rank = {v: i for i, v in enumerate(word, 1)}
    sing = {r: False for r in roots}
    for basin in first_merge(word, masks):
        if basin.bit_count() == 1:
            sing[basin.bit_length()-1] = True
    found = []
    for rescued, erased in ((p, q), (q, p)):
        for v in (erased, *neighbors[erased]):
            assert v != rescued and rescued not in neighbors[v]
            max_neighbor = min(neighbors[v], key=rank.get)
            erased_neighbor_peaks = [
                u for u in neighbors[v] if all(rank[u] < rank[w] for w in neighbors[u] if w != v)
            ]
            t_rank = min((rank[u] for u in erased_neighbor_peaks), default=len(word)+1)
            gap = t_rank-rank[max_neighbor]
            expected = sing[erased] if v == erased else rank[max_neighbor] < rank[erased]
            assert (gap > 0) == expected
            if v in neighbors[erased] and gap > 0:
                assert erased == q and t_rank == rank[erased]
            if gap > 0:
                found.append(gap)
    lens = any(distance(v, q, points) == 2 for v in word[:rank[q]-1])
    assert bool(found) == (any(sing.values()) or lens)
    return d, bool(found)


def classify():
    for n in (3, 4, 5):
        points, _, neighbors, masks = grid(n)
        rng = Random(20261118+n)
        words = permutations(range(n*n)) if n == 3 else (rng.sample(range(n*n), n*n) for _ in range(25000))
        results = Counter()
        for word in words:
            ans = check_word(word, points, neighbors, masks)
            if ans:
                results[ans] += 1
        print(f"n={n} exact per-word geometry checks (distance, eligible): {sorted(results.items())}")


def audit_path(d):
    n = d+5
    points, ids, _, masks = grid(n)
    p, q = ids[(1,2)], ids[(1+d,2)]
    interior = tuple(ids[(x,2)] for x in range(2,d+1))
    strip = tuple(ids[(x,y)] for x in range(2,d+1) for y in (1,3))
    outside = set(range(n*n)) - {p,q,*interior,*strip}
    root = ids[(0,2)]
    bfs = deque([root])
    seen = {root}
    order = [root]
    while bfs:
        v = bfs.popleft()
        more = [u for u in outside-seen if masks[v] & (1<<u)]
        for u in sorted(more):
            seen.add(u)
            order.append(u)
            bfs.append(u)
    assert seen == outside
    base = {u: F(1,2) + F(len(order)-i, 7*len(order)) for i,u in enumerate(order)}
    base.update({u: F(1,4) + F(i+1, 20*len(strip)) for i,u in enumerate(strip)})
    N = n*n
    a, b = 1-F(1,3*N), 1-F(3,2*N)
    base[p], base[q] = a, b
    for blocks in product((0,1,2), repeat=len(interior)):
        values = base.copy()
        for i,(u,k) in enumerate(zip(interior, blocks), 1):
            shelves = [(F(0),F(1,4)), (F(1,3),b), (b,a)]
            l,h = shelves[k]
            values[u] = l + (h-l)*F(i,d)
        if not all(k == 2 for k in blocks):
            assert len(peaks(sorted(values, key=values.get, reverse=True), masks)) >= 2
        if not any(blocks):
            assert peaks(sorted(values, key=values.get, reverse=True), masks) == [p,q]
    increasing = base | {u: a - (a-b)*F(i,d) for i,u in enumerate(interior,1)}
    assert peaks(sorted(increasing, key=increasing.get, reverse=True), masks) == [p]
    return 3**(d-1)


if __name__ == "__main__":
    classify()
    cases = sum(audit_path(d) for d in range(3,9))
    print(f"fixed-path fibre test inputs: {cases}")
    for d in range(3,101):
        assert sum(4*j for j in (d-1,d,d+1))*2 == 24*d
    print("near-anchor shell multiplicity formula: d=3..100")
