#!/usr/bin/env python3
"""Test the structured block-transpose relabelling for R^a D^b paths."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = ('r', 'r')
HERE = Path(__file__).resolve().parent


def grid_edges(n):
    for x in range(n):
        for y in range(n):
            if x + 1 < n:
                yield (x, y), (x + 1, y)
            if y + 1 < n:
                yield (x, y), (x, y + 1)


def quotient(n, S):
    verts = {ROOT} | {(x, y) for x in range(n) for y in range(n) if (x, y) not in S}
    order = [ROOT] + sorted(verts - {ROOT})
    pos = {v: i for i, v in enumerate(order)}
    adj = [0] * len(order)
    for u, v in grid_edges(n):
        u = ROOT if u in S else u; v = ROOT if v in S else v
        if u != v:
            i, j = pos[u], pos[v]
            adj[i] |= 1 << j; adj[j] |= 1 << i
    return order, tuple(adj)


def separated(G, blockers, target):
    seen = todo = 1
    allowed = ((1 << len(G)) - 1) & ~blockers
    while todo:
        bit = todo & -todo; todo -= bit
        i = bit.bit_length() - 1
        new = G[i] & allowed & ~seen
        seen |= new; todo |= new
    return not (seen >> target & 1)


def avoiding_path(G, blockers, target):
    allowed = ((1 << len(G)) - 1) & ~blockers
    todo, parent = [0], {0: None}
    while todo:
        i = todo.pop(0)
        if i == target:
            out = []
            while i is not None:
                out.append(i); i = parent[i]
            return out[::-1]
        candidates = G[i] & allowed
        for j in range(len(G)):
            if candidates >> j & 1 and j not in parent:
                parent[j] = i; todo.append(j)
    return None


def mapped_graph(G, perm):
    ans = [0] * len(G)
    for i in range(len(G)):
        for j in range(i + 1, len(G)):
            if G[i] >> j & 1:
                a, b = perm[i], perm[j]
                ans[a] |= 1 << b; ans[b] |= 1 << a
    return tuple(ans)


def path(a, b):
    return {(x, 0) for x in range(a + 1)} | {(a, y) for y in range(1, b + 1)}


def one(n, a):
    b = n - 1 - a
    old_order, G = quotient(n, {(x, 0) for x in range(n)})
    new_order, H = quotient(n, path(a, b))
    newpos = {v: i for i, v in enumerate(new_order)}
    perm = [0]
    images = []
    for v in old_order[1:]:
        x, y = v
        image = (a + y, x - a) if a <= x <= a + b and 1 <= y <= b else v
        assert image in newpos
        images.append(image); perm.append(newpos[image])
    assert len(set(perm)) == len(perm)
    MG = mapped_graph(G, perm)
    bad = []
    for v in range(1, len(H)):
        if H[0] >> v & 1:
            continue
        if not separated(MG, H[v], v):
            witness = avoiding_path(MG, H[v], v)
            assert witness is not None
            bad.append({'vertex': new_order[v], 'path': [new_order[i] for i in witness]})
    return {'n': n, 'a': a, 'b': b, 'word': 'R' * a + 'D' * b,
            'block_transpose_is_bijection': True,
            'cut_criterion_passed': not bad,
            'failed_target_vertices': [{'vertex': list(x['vertex']),
                                        'avoiding_path': ['root' if v == ROOT else list(v) for v in x['path']]}
                                       for x in bad]}


def main():
    rows = [one(n, a) for n in range(3, 10) for a in range(1, n - 1)]
    out = {'kind': 'one-turn-block-transpose-cut-audit',
           'scope': 'Exact bitset test of one explicit relabelling only. A failing row refutes this proposed relabelling, not all relabellings.',
           'rows': rows,
           'passed': sum(x['cut_criterion_passed'] for x in rows),
           'total': len(rows)}
    text = json.dumps(out, indent=2) + '\n'
    HERE.joinpath('one-turn-block-transpose-audit.json').write_text(text)
    print(text, end='')
    print('output_sha256=' + hashlib.sha256(text.encode()).hexdigest())


if __name__ == '__main__':
    main()
