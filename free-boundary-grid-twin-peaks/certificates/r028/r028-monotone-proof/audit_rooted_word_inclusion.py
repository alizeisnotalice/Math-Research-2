#!/usr/bin/env python3
"""Exhaustive finite audit of the rooted-word separation criterion."""
from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent


def masks_from_code(n: int, code: int) -> tuple[int, ...]:
    a = [0] * n
    bit = 0
    for i in range(n):
        for j in range(i + 1, n):
            if code >> bit & 1:
                a[i] |= 1 << j
                a[j] |= 1 << i
            bit += 1
    return tuple(a)


def connected(a: tuple[int, ...]) -> bool:
    seen, todo = 1, 1
    while todo:
        bit = todo & -todo
        todo -= bit
        v = bit.bit_length() - 1
        new = a[v] & ~seen
        seen |= new
        todo |= new
    return seen == (1 << len(a)) - 1


def rooted_word_mask(a: tuple[int, ...]) -> int:
    n, out = len(a), 0
    for index, tail in enumerate(itertools.permutations(range(1, n))):
        seen, ok = 1, True
        for v in tail:
            if not (a[v] & seen):
                ok = False; break
            seen |= 1 << v
        if ok:
            out |= 1 << index
    return out


def separated(a: tuple[int, ...], blockers: int, target: int) -> bool:
    if blockers & 1:
        return True
    seen, todo = 1, 1
    allowed = ((1 << len(a)) - 1) & ~blockers
    while todo:
        bit = todo & -todo
        todo -= bit
        v = bit.bit_length() - 1
        new = a[v] & allowed & ~seen
        seen |= new
        todo |= new
    return not (seen >> target & 1)


def criterion(G: tuple[int, ...], H: tuple[int, ...]) -> bool:
    for v in range(1, len(G)):
        if H[0] >> v & 1:
            continue
        if not separated(G, H[v], v):
            return False
    return True


def audit(n: int) -> dict[str, object]:
    m = n * (n - 1) // 2
    graphs = [masks_from_code(n, c) for c in range(1 << m)]
    graphs = [a for a in graphs if connected(a)]
    word_masks = [rooted_word_mask(a) for a in graphs]
    checks = 0
    for i, G in enumerate(graphs):
        for j, H in enumerate(graphs):
            direct = word_masks[i] & ~word_masks[j] == 0
            cut = criterion(G, H)
            assert direct == cut
            checks += 1
    return {'vertices': n, 'connected_rooted_graphs': len(graphs),
            'ordered_graph_pairs_checked': checks, 'all_passed': True}


def main() -> None:
    out = {'kind': 'exhaustive-rooted-word-inclusion-cut-criterion-audit',
           'scope': 'All connected simple graphs on a fixed rooted labelled vertex set, with root 0. Direct word containment is compared to the separation criterion.',
           'runs': [audit(3), audit(4), audit(5)]}
    text = json.dumps(out, indent=2) + '\n'
    HERE.joinpath('rooted-word-inclusion-audit.json').write_text(text)
    print(text, end='')
    print('output_sha256=' + hashlib.sha256(text.encode()).hexdigest())


if __name__ == '__main__':
    main()
