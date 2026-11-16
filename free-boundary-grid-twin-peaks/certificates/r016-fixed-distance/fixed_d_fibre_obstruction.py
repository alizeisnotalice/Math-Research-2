#!/usr/bin/env python3
"""Exact checks for a fixed-distance-three large canonical lifting fibre.

For n>=5, roots q=(0,0), p=(3,0) have their unique geodesic through
w1=(1,0), w2=(2,0).  A connected H-prefix avoids p and its H-neighbours,
then p is born.  The two removed path points can be inserted in arbitrary
order and arbitrary positions after p, all mapping under path promotion to
the same one-peak word.
"""
from itertools import combinations, permutations
from pathlib import Path
import json
import platform
import sys


def neighbours(v, n):
    x, y = v
    return {(X, Y) for X, Y in ((x-1,y),(x+1,y),(x,y-1),(x,y+1))
            if 0 <= X < n and 0 <= Y < n}


def peaks(word, n):
    pos = {v: i for i, v in enumerate(word)}  # descending-label ranks
    return tuple(v for v in word
                 if all(pos[v] < pos[x] for x in neighbours(v, n) if x in pos))


def grow_order(allowed, seed, n, stop=None):
    order = list(seed)
    seen = set(seed)
    assert len(order) == len(seen) and seen <= allowed
    while len(order) < (len(allowed) if stop is None else stop):
        choices = sorted(v for v in allowed-seen if neighbours(v, n) & seen)
        assert choices
        order.append(choices[0])
        seen.add(choices[0])
    return order


def connected(vertices, n):
    todo = [next(iter(vertices))]
    seen = set(todo)
    while todo:
        v = todo.pop()
        for x in neighbours(v, n) & vertices:
            if x not in seen:
                seen.add(x)
                todo.append(x)
    return seen == vertices


def instance(n):
    assert n >= 5
    V = {(x, y) for x in range(n) for y in range(n)}
    q, p = (0, 0), (3, 0)
    path = ((1, 0), (2, 0))
    anchors = ((1, 1), (2, 1))
    H = V-set(path)
    forbidden = {p} | (neighbours(p, n) & H)
    K = H-forbidden
    assert connected(K, n)
    seed = (q, (0, 1), *anchors)
    r = (n*n)//2
    assert r >= len(seed) and not (neighbours(p, n) & set(seed))
    pre = grow_order(K, seed, n, stop=r)
    # Finish H after p; all future H points attach to the current H set.
    Hword = grow_order(H, tuple(pre)+(p,), n)
    assert Hword[:r] == pre and Hword[r] == p
    assert peaks(Hword, n) == (q, p)
    tail = Hword[r+1:]
    h = len(tail)
    expected = (h+1)*(h+2)  # choose two slots and their order
    common = None
    checked = 0
    for slots in combinations(range(h+2), 2):
        for pp in permutations(path):
            it_tail = iter(tail); it_path = iter(pp)
            merged = tuple(next(it_path) if i in slots else next(it_tail)
                           for i in range(h+2))
            source = tuple(pre)+(p,)+merged
            assert peaks(source, n) == (q, p)
            output = (q,)+path+tuple(Hword[1:])
            assert peaks(output, n) == (q,)
            assert output[:3] == (q,)+path
            if common is None:
                common = output
            assert output == common
            checked += 1
    assert checked == expected
    # q to p is horizontal, hence has a unique graph geodesic.
    assert abs(q[0]-p[0]) + abs(q[1]-p[1]) == 3
    return {'n': n, 'N': n*n, 'roots': {'q': q, 'p': p}, 'path_internal': path,
            'prefix_size_before_p': r, 'tail_size_after_p': h,
            'K_connected': True,
            'same_output_true_fibre_lower_bound': expected,
            'checked_sources': checked, 'all_exactly_two_peaks': True,
            'common_output_is_one_peak': True,
            'unique_root_geodesic': True}


out = {'status': 'passed', 'python': sys.version, 'platform': platform.platform(),
       'arithmetic': 'finite deterministic enumeration; no random seed',
       'rows': [instance(n) for n in (5, 6, 7, 8)]}
Path('fixed-d-fibre-obstruction.json').write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps(out, indent=2))
