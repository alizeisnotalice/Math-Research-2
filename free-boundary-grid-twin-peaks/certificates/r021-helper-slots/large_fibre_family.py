#!/usr/bin/env python3
"""Full-grid Theta(N) fibre family for the single-slot stable-helper map."""
import json


def family(n):
    assert n >= 5
    p, r, q, z = (1, 2), (2, 1), (4, 2), (2, 2)
    V = [(x, y) for y in range(n) for x in range(n)]
    neigh = {v: [(xx, yy) for xx, yy in
                 ((v[0]-1, v[1]), (v[0]+1, v[1]),
                  (v[0], v[1]-1), (v[0], v[1]+1))
             if 0 <= xx < n and 0 <= yy < n]
             for v in V}
    roots = [p, r, q]
    seen, queue, tail = set(roots), roots.copy(), []
    for v in queue:
        for u in neigh[v]:
            if u != z and u not in seen:
                seen.add(u); queue.append(u); tail.append(u)
    assert seen == set(V) - {z}
    base = roots + tail  # a legal three-root word after z is deleted
    w = (3, 2)
    assert p in neigh[z] and r in neigh[z] and q not in neigh[z]
    common = set(V)
    for v in roots:
        common &= set(neigh[v])
    assert not common
    closed = {v: set(neigh[v]) | {v} for v in V}
    ordered_cover = [(a, b) for a in V for b in neigh[a]
                     if set(roots) <= (closed[a] | closed[b])]
    assert (z, w) in ordered_cover and (w, z) in ordered_cover
    assert len(ordered_cover) == 2

    def peakset(word):
        earlier, ans = set(), []
        for v in word:
            if not (set(neigh[v]) & earlier):
                ans.append(v)
            earlier.add(v)
        return ans

    outputs = set()
    # Slots strictly after q in base: N-3 choices.
    for slot in range(3, len(base) + 1):
        source = base.copy(); source.insert(slot, z)
        assert peakset(source) == roots
        target = source.copy(); target.remove(z); target.insert(1, z)
        assert peakset(target) == [p, q]
        outputs.add(tuple(target))
    assert len(outputs) == 1
    return {'n': n, 'source_words_with_same_marked_output': n*n-3,
            'source_peakset': roots, 'target_peakset': [p, q],
            'common_neighbour_count': len(common),
            'ordered_domino_cover_count': len(ordered_cover),
            'marked_vertices': {'p': p, 'z': z, 'q': q, 'r': r}}


print(json.dumps({'status': 'passed',
                  'families': [family(n) for n in range(5, 13)]}, indent=2))
