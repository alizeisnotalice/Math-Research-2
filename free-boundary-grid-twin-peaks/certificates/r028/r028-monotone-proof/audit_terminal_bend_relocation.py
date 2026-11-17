#!/usr/bin/env python3
"""Exact audit of the top-row to terminal-bend rooted-graph transformation."""
from __future__ import annotations

import hashlib
import json
from functools import lru_cache
from pathlib import Path


ROOT = ('r', 'r')
HERE = Path(__file__).resolve().parent


def grid_edges(n: int):
    ans = set()
    for x in range(n):
        for y in range(n):
            if x + 1 < n:
                ans.add(frozenset(((x, y), (x + 1, y))))
            if y + 1 < n:
                ans.add(frozenset(((x, y), (x, y + 1))))
    return ans


def quotient(n: int, S: set[tuple[int, int]]):
    verts = {ROOT} | {(x, y) for x in range(n) for y in range(n) if (x, y) not in S}
    edges = set()
    for e in grid_edges(n):
        u, v = tuple(e)
        u = ROOT if u in S else u
        v = ROOT if v in S else v
        if u != v:
            edges.add(frozenset((u, v)))
    return verts, edges


def count_orders(verts, edges):
    order = [ROOT] + sorted(v for v in verts if v != ROOT)
    pos = {v: i for i, v in enumerate(order)}
    adj = [0] * len(order)
    for e in edges:
        u, v = tuple(e)
        i, j = pos[u], pos[v]
        adj[i] |= 1 << j
        adj[j] |= 1 << i
    full = (1 << len(order)) - 1
    @lru_cache(None)
    def h(mask):
        if mask == full:
            return 1
        total = 0
        for i in range(1, len(order)):
            if not (mask >> i & 1) and adj[i] & mask:
                total += h(mask | (1 << i))
        return total
    return h(1), order, adj


def valid(word, adj):
    seen = 1
    for i in word[1:]:
        if not (adj[i] & seen):
            return False
        seen |= 1 << i
    return True


def stream_verify_injection(old_order, old_adj, new_order, new_adj, phi):
    """Visit every old word without retaining any word list in memory."""
    full = (1 << len(old_order)) - 1
    new_pos = {v: i for i, v in enumerate(new_order)}
    images = [new_pos[phi[v]] for v in old_order]
    checked = 0
    def rec(old_mask, new_mask):
        nonlocal checked
        if old_mask == full:
            checked += 1
            return
        for i in range(1, len(old_order)):
            if not (old_mask >> i & 1) and old_adj[i] & old_mask:
                j = images[i]
                assert new_adj[j] & new_mask
                rec(old_mask | (1 << i), new_mask | (1 << j))
    rec(1, 1)
    return checked


def root_then_z_completion_exists(order, adj, z):
    zindex = order.index(z)
    assert adj[zindex] & 1
    seen = (1 << 0) | (1 << zindex)
    full = (1 << len(order)) - 1
    while seen != full:
        options = [i for i in range(1, len(order)) if not (seen >> i & 1) and adj[i] & seen]
        assert options
        seen |= 1 << options[0]
    return True


def one(n: int):
    T = {(x, 0) for x in range(n)}
    S = {(x, 0) for x in range(n - 1)} | {(n - 2, 1)}
    old_v, old_e = quotient(n, T)
    new_v, new_e = quotient(n, S)
    oldx, newy = (n - 2, 1), (n - 1, 0)
    a, b, z = (n - 3, 1), (n - 1, 1), (n - 2, 2)
    phi = {ROOT: ROOT, oldx: newy}
    phi.update({v: v for v in old_v if v not in (ROOT, oldx)})
    assert set(phi.values()) == new_v
    mapped = {frozenset(phi[v] for v in e) for e in old_e}
    drop_a, drop_z, add_z = frozenset((newy, a)), frozenset((newy, z)), frozenset((ROOT, z))
    assert mapped - new_e == {drop_a, drop_z}
    assert new_e - mapped == {add_z}
    assert frozenset((ROOT, newy)) in mapped and frozenset((ROOT, a)) in mapped
    direct = {'n': n,
              'mapped_edge_removals': [sorted(map(str, e)) for e in (drop_a, drop_z)],
              'root_edge_addition': sorted(map(str, add_z))}
    if n in (3, 4):
        old_count, old_order, old_adj = count_orders(old_v, old_e)
        new_count, new_order, new_adj = count_orders(new_v, new_e)
        assert old_count < new_count
        direct.update({'top_h': str(old_count), 'terminal_bend_h': str(new_count),
                       'strict': old_count < new_count})
        # Convert every old basic word through phi and verify it is a new word.
        checked = stream_verify_injection(old_order, old_adj, new_order, new_adj, phi)
        assert checked == old_count
        assert root_then_z_completion_exists(new_order, new_adj, z)
        direct['all_old_words_injected_into_new'] = True
        direct['old_words_stream_checked'] = checked
        direct['new_word_starting_root_then_z_exists'] = True
    return direct


def main():
    out = {'kind': 'terminal-bend-root-relocation-audit',
           'scope': 'Exact rooted graph transformation for T=R^(n-1) and S=R^(n-2)D. Full rooted-word injection is exhaustively checked only n=3,4; the general proof is separate.',
           'runs': [one(n) for n in range(3, 8)]}
    text = json.dumps(out, indent=2) + '\n'
    HERE.joinpath('terminal-bend-relocation-audit.json').write_text(text)
    print(text, end='')
    print('output_sha256=' + hashlib.sha256(text.encode()).hexdigest())


if __name__ == '__main__':
    main()
