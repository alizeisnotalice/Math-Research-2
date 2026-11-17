#!/usr/bin/env python3
"""Exact audit of the epsilon-free cylinders in the rank-two graph case."""
import hashlib
import itertools
import json
import platform
import time
from pathlib import Path
from build_graph_cases import basis, fwht


def parity(x):
    return x.bit_count() & 1


def span(bs):
    s = {0}
    for b in bs:
        s |= {x ^ b for x in tuple(s)}
    return s


def alt_rank(form, n=6):
    rows = [0] * n
    for j, (a, b) in enumerate(itertools.combinations(range(n), 2)):
        if (form >> j) & 1:
            rows[a] ^= 1 << b
            rows[b] ^= 1 << a
    return len(basis(rows))


def main():
    t0 = time.monotonic()
    # Physical bits y0,y1,y2,t; parameter bits c0,l0,l1,l2,lt,hQ.
    quad = lambda x: ((x >> 0) & 1) & ((x >> 1) & 1)
    base = lambda x: (((x >> 0) & 1) & ((x >> 2) & 1)) ^ (((x >> 3) & 1) & ((x >> 1) & 1)) ^ quad(x)
    blocks = [[x for x in range(16) if (((x >> 3) & 1) ^ quad(x)) == c] for c in (0, 1)]
    generators = [(1 << 16) - 1] + [sum(((x >> j) & 1) << x for x in range(16)) for j in range(4)] + [sum(quad(x) << x for x in range(16))]
    assert len(basis(generators)) == 6
    cylinders = {(c, u): set() for c in (0, 1) for u in (0, 1)}
    for p in range(64):
        phase = sum(base(x) << x for x in range(16))
        for j, g in enumerate(generators):
            if (p >> j) & 1:
                phase ^= g
        assert all(abs(w) == 4 for w in fwht([1 - 2 * ((phase >> x) & 1) for x in range(16)]))
        for c, block in enumerate(blocks):
            support = {i for i, w in enumerate(fwht([(1 - 2 * ((phase >> x) & 1)) if x in block else 0 for x in range(16)])) if w}
            predicted = set()
            for u in (0, 1):
                for xi in range(16):
                    eq1 = ((xi >> 2) & 1) ^ ((p >> 3) & 1)
                    eq2 = ((xi >> 1) & 1) ^ ((p >> 2) & 1) ^ (u & (((p >> 4) ^ (p >> 5) ^ (xi >> 3)) & 1))
                    if eq1 == u and eq2 == c:
                        cylinders[c, u].add(p | (xi << 6))
                        predicted.add(xi)
            assert support == predicted
    covered = set()
    for s in cylinders.values():
        assert not covered & s
        covered |= s
        a = min(s)
        directions = [x ^ a for x in s]
        assert len(basis(directions)) == 8 and len(s) == 256
    assert covered == set(range(1024))
    normals = [(1 << 3) ^ (1 << 8), (1 << 2) ^ (1 << 7),
               (1 << 2) ^ (1 << 7) ^ (1 << 4) ^ (1 << 5) ^ (1 << 9)]
    assert len(basis(normals)) == len(basis([n & 63 for n in normals])) == 3
    kernel = [p for p in range(64) if all(parity(p & n) == 0 for n in normals)]
    assert len(kernel) == 8
    kranks = []
    for block in blocks:
        patterns = []
        for p in kernel:
            g = 0
            for j, b in enumerate(generators):
                if (p >> j) & 1:
                    g ^= b
            patterns.append(sum(((g >> x) & 1) << j for j, x in enumerate(block)))
        kranks.append(len(basis(patterns)))
    assert kranks == [2, 2]
    edges = list(itertools.combinations(range(6), 2))
    forms = [f for f in range(1 << 15) if alt_rank(f) == 2]
    compatible = {}
    for rank_q in (2, 4, 6):
        qform = sum(1 << edges.index((j, j + 1)) for j in range(0, rank_q, 2))
        compatible[str(rank_q)] = sum(alt_rank(f ^ qform) == 2 for f in forms)
    assert compatible == {'2': 90, '4': 20, '6': 0}
    root = Path(__file__).parent
    result = {'schema': 'a28-r018-rank2-graph-cylinders-v1', 'phase_count': 64,
              'joint_points_checked': 1024, 'cylinder_count': 4, 'codimension': 2,
              'normal_rank': 3, 'projected_rank': 3, 'kernel_dimension': 3,
              'kernel_block_ranks': kranks, 'six_dimensional_rank2_forms': len(forms),
              'rank2_compatible_with_fixed_Q_by_rank': compatible,
              'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'helper_sha256': hashlib.sha256((root / 'build_graph_cases.py').read_bytes()).hexdigest(),
              'python': platform.python_version(), 'elapsed_seconds': time.monotonic() - t0,
              'arithmetic': 'exact GF(2), integer Walsh', 'failed_assertions': 0}
    (root / 'rank2-graph-cylinders.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
