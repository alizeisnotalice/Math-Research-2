#!/usr/bin/env python3
"""Verify the turn loss formula for contracted monotone-path root graphs."""
from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent


def path_cells(n: int, word: str) -> set[tuple[int, int]]:
    x = y = 0
    out = {(x, y)}
    for c in word:
        x, y = (x + 1, y) if c == 'R' else (x, y + 1)
        out.add((x, y))
    return out


def grid_edges(n: int) -> set[tuple[tuple[int, int], tuple[int, int]]]:
    ans = set()
    for x in range(n):
        for y in range(n):
            if x + 1 < n:
                ans.add(((x, y), (x + 1, y)))
            if y + 1 < n:
                ans.add(((x, y), (x, y + 1)))
    return ans


def contracted_edge_count(n: int, S: set[tuple[int, int]]) -> int:
    # Contract S to a root and delete loops / merge repeated root edges.
    root = ('root', 'root')
    edges = set()
    for u, v in grid_edges(n):
        uu = root if u in S else u
        vv = root if v in S else v
        if uu != vv:
            edges.add(tuple(sorted((uu, vv), key=repr)))
    return len(edges)


def main() -> None:
    runs = []
    examples = []
    for n in range(3, 8):
        E = 2 * n * (n - 1)
        rows = []
        for bits in itertools.product('RD', repeat=n - 1):
            word = ''.join(bits)
            turns = sum(word[i] != word[i + 1] for i in range(len(word) - 1))
            actual = contracted_edge_count(n, path_cells(n, word))
            expected = E - (n - 1) - turns
            assert actual == expected
            rows.append({'word': word, 'turns': turns, 'contracted_simple_edges': actual})
        runs.append({'n': n, 'grid_edges': E, 'path_count': len(rows),
                     'all_turn_formula_passed': True})
        if n == 5:
            examples = [x for x in rows if x['word'] in ('RRRR', 'RDRD', 'DDRD')]
    out = {
        'kind': 'exact-contracted-root-graph-turn-loss-audit',
        'formula': '|E(G_n/S_w)| = 2n(n-1) - (n-1) - turns(w)',
        'scope': 'Pure graph calculation for every corner monotone path n=3..7; it checks no completion counts.',
        'runs': runs, 'n5_examples': examples,
    }
    text = json.dumps(out, indent=2) + '\n'
    HERE.joinpath('contraction-turn-edges.json').write_text(text)
    print(text, end='')
    print('output_sha256=' + hashlib.sha256(text.encode()).hexdigest())


if __name__ == '__main__':
    main()
