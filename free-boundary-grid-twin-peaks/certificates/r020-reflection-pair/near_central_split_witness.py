#!/usr/bin/env python3
"""Construct and verify central-split two-peak words at distance two."""
import json


def witness(n):
    m = n // 2
    p = (m - 1, 0)
    q = (m, 1)
    vertices = [(x, y) for y in range(n) for x in range(n)]
    neighbours = {v: [(xx, yy) for xx, yy in
                      ((v[0] - 1, v[1]), (v[0] + 1, v[1]),
                       (v[0], v[1] - 1), (v[0], v[1] + 1))
                      if 0 <= xx < n and 0 <= yy < n]
                  for v in vertices}
    seen, queue, tail = {q}, [q], []
    for v in queue:
        for u in neighbours[v]:
            if u != p and u not in seen:
                seen.add(u); queue.append(u); tail.append(u)
    assert seen == set(vertices) - {p}
    word = [q, p] + tail
    assert len(word) == n * n and len(set(word)) == n * n
    earlier, peaks = set(), []
    for v in word:
        if not (set(neighbours[v]) & earlier): peaks.append(v)
        earlier.add(v)
    assert peaks == [q, p]
    assert abs(p[0] - q[0]) + abs(p[1] - q[1]) == 2
    assert p[0] < m <= q[0]
    return {'n': n, 'p': p, 'q': q, 'peaks': peaks,
            'distance': 2, 'word_high_to_low': word}


print(json.dumps({'status': 'passed',
                  'witnesses': [witness(n) for n in range(2, 21, 2)]}, indent=2))
