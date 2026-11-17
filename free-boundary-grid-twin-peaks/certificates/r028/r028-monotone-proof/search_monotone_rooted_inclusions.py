#!/usr/bin/env python3
"""Bounded deterministic search for root-preserving inclusion relabellings.

This is a witness finder only.  A failed bounded search is explicitly not a
nonexistence certificate.  Every returned permutation is certified by the
separation criterion from rooted_word_inclusion.md.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import random
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
    vertices = {ROOT} | {(x, y) for x in range(n) for y in range(n) if (x, y) not in S}
    edges = set()
    for u, v in grid_edges(n):
        u = ROOT if u in S else u
        v = ROOT if v in S else v
        if u != v:
            edges.add(tuple(sorted((u, v), key=repr)))
    order = [ROOT] + sorted(vertices - {ROOT})
    pos = {v: i for i, v in enumerate(order)}
    adj = [0] * len(order)
    for u, v in edges:
        i, j = pos[u], pos[v]
        adj[i] |= 1 << j; adj[j] |= 1 << i
    return order, tuple(adj)


def separated(G, blockers, target):
    seen, todo = 1, 1
    allowed = ((1 << len(G)) - 1) & ~blockers
    while todo:
        bit = todo & -todo; todo -= bit
        i = bit.bit_length() - 1
        new = G[i] & allowed & ~seen
        seen |= new; todo |= new
    return not (seen >> target & 1)


def mapped_graph(G, perm):
    out = [0] * len(G)
    for i in range(len(G)):
        for j in range(i + 1, len(G)):
            if G[i] >> j & 1:
                a, b = perm[i], perm[j]
                out[a] |= 1 << b; out[b] |= 1 << a
    return tuple(out)


def failures(Gmapped, H):
    bad = []
    for v in range(1, len(H)):
        if H[0] >> v & 1:
            continue
        if not separated(Gmapped, H[v], v):
            bad.append(v)
    return bad


def path_set(n, word):
    x = y = 0; S = {(x, y)}
    for c in word:
        x, y = (x + 1, y) if c == 'R' else (x, y + 1)
        S.add((x, y))
    return S


def common_start(source_labels, new_order):
    newpos = {v: i for i, v in enumerate(new_order)}
    perm = [0] * len(source_labels)
    missing_old, missing_new = [], []
    for i, v in enumerate(source_labels[1:], 1):
        if v in newpos:
            perm[i] = newpos[v]
        else:
            missing_old.append(i)
    used = set(perm[1:])
    missing_new = [i for i in range(1, len(new_order)) if i not in used]
    for i, j in zip(missing_old, missing_new):
        perm[i] = j
    return perm


def find_witness(G, H, initials, seed, restarts=24, steps=1200):
    def score(p): return len(failures(mapped_graph(G, p), H))
    best, bestscore = initials[0][1][:], score(initials[0][1])
    for label, initial in initials:
        initial_score = score(initial)
        if initial_score < bestscore:
            best, bestscore = initial[:], initial_score
        if initial_score == 0:
            return initial, {'method': label, 'score': 0}
    rng = random.Random(seed)
    for restart in range(restarts):
        p = initials[restart % len(initials)][1][:]
        if restart:
            tail = p[1:]; rng.shuffle(tail); p[1:] = tail
        current = score(p)
        for step in range(steps):
            i, j = rng.sample(range(1, len(p)), 2)
            p[i], p[j] = p[j], p[i]
            trial = score(p)
            # Greedy descent with a small deterministic temperature for plateaux.
            if trial <= current or (trial == current + 1 and (step + restart) % 97 == 0):
                current = trial
            else:
                p[i], p[j] = p[j], p[i]
            if current < bestscore:
                best, bestscore = p[:], current
            if current == 0:
                return p, {'method': 'bounded seeded swap search', 'restart': restart, 'step': step, 'score': 0}
    return None, {'method': 'bounded seeded swap search', 'best_score': bestscore,
                  'restarts': restarts, 'steps_per_restart': steps}


def one(n, word):
    top = {(x, 0) for x in range(n)}
    old_order, G = quotient(n, top)
    new_order, H = quotient(n, path_set(n, word))
    initial = common_start(old_order, new_order)
    transposed = [ROOT] + [(y, x) for x, y in old_order[1:]]
    initials = [('common-coordinate start', initial),
                ('transpose then common-coordinate start', common_start(transposed, new_order))]
    seed = 20261117 + 1000 * n + sum((i + 1) * (c == 'D') for i, c in enumerate(word))
    perm, detail = find_witness(G, H, initials, seed)
    result = {'n': n, 'word': word, 'vertices_after_contraction': len(G), 'search': detail}
    if perm is not None:
        GM = mapped_graph(G, perm)
        assert not failures(GM, H)
        result['status'] = 'witness_found'
        result['mapping'] = [{'from': list(old_order[i]) if old_order[i] != ROOT else 'root',
                              'to': list(new_order[perm[i]]) if new_order[perm[i]] != ROOT else 'root'}
                             for i in range(len(perm))]
    else:
        result['status'] = 'no_witness_within_bounded_search'
    return result


def main():
    runs = []
    for n in (4, 5):
        for bits in itertools.product('RD', repeat=n - 1):
            runs.append(one(n, ''.join(bits)))
    out = {'kind': 'bounded-rooted-inclusion-relabel-search',
           'scope': 'Witness search only; no-witness rows do not certify nonexistence.',
           'found': sum(x['status'] == 'witness_found' for x in runs),
           'total': len(runs), 'runs': runs}
    text = json.dumps(out, indent=2) + '\n'
    HERE.joinpath('monotone-rooted-inclusion-search.json').write_text(text)
    print(json.dumps({'found': out['found'], 'total': out['total'],
                      'statuses': {str(n): [x['status'] for x in runs if x['n'] == n] for n in (4, 5)}}, indent=2))
    print('output_sha256=' + hashlib.sha256(text.encode()).hexdigest())


if __name__ == '__main__':
    main()
