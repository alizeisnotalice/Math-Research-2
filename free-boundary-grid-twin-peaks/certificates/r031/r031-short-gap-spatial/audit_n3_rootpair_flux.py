#!/usr/bin/env python3
"""Full 9! audit of the directed and symmetrized h2-2h3 root-pair kernel."""

from collections import defaultdict
from fractions import Fraction
from itertools import permutations
import json
from pathlib import Path


N = 3
V = tuple(range(N * N))
ADJ = {
    v: {yy * N + xx for xx, yy in ((v % N + 1, v // N), (v % N - 1, v // N),
                                     (v % N, v // N + 1), (v % N, v // N - 1))
        if 0 <= xx < N and 0 <= yy < N}
    for v in V
}


def geometry(u, v):
    dx, dy = abs(u % N - v % N), abs(u // N - v // N)
    return dx + dy, dx * dx + dy * dy


def main():
    h = {2: defaultdict(Fraction), 3: defaultdict(Fraction)}
    double_peak = 0
    for word in permutations(V):
        pos = [0] * len(V)
        for i, x in enumerate(word):
            pos[x] = i
        peaks = [x for x in V if all(pos[x] < pos[y] for y in ADJ[x])]
        if len(peaks) != 2:
            continue
        double_peak += 1
        v = word[0]
        p = peaks[1] if peaks[0] == v else peaks[0]
        z = [x for x in word if x != v]
        z.remove(p)
        z.insert(0, p)
        L = next(i for i, x in enumerate(z, 1) if x in ADJ[p]) - 1
        if L in h:
            h[L][v, p] += Fraction(1, L)
    assert double_peak == 143808
    assert all(x.denominator == 1 for table in h.values() for x in table.values())
    h_int = {L: {f'{v},{p}': int(x) for (v, p), x in table.items()} for L, table in h.items()}
    c = defaultdict(Fraction)
    for pair in set(h[2]) | set(h[3]):
        c[pair] = h[2][pair] - 2 * h[3][pair]
    assert all(x.denominator == 1 for x in c.values())
    c = {pair: int(value) for pair, value in c.items() if value}
    # First audit against the archived distance/e2 aggregate.
    aggregate = defaultdict(int)
    for (v, p), value in c.items():
        aggregate[geometry(v, p)] += value
    reference = json.loads(Path(__file__).with_name('joint_h2_h3_audit.json').read_text())
    r3 = reference['results'][0]
    expected = {(cell['d'], cell['e2']): int(cell['h2_minus_2h3'])
                for cell in r3['nonzero_joint_cells']}
    assert {k: v for k, v in aggregate.items() if v} == expected

    # Symmetrize the directed deletion orientation.  It is the natural kernel
    # against a symmetric root-pair observable F({u,v}).
    sym = {}
    for u in V:
        for v in V:
            if u < v:
                z = c.get((u, v), 0) + c.get((v, u), 0)
                if z:
                    sym[u, v] = z
    row_sums = {u: sum(z for (a, b), z in sym.items() if u in (a, b)) for u in V}
    out = {
        'status': 'passed',
        'scope': 'complete n=3 finite root-pair audit only',
        'directed_h2': h_int[2],
        'directed_h3': h_int[3],
        'directed_delta_nonzero': {f'{v},{p}': z for (v, p), z in sorted(c.items())},
        'symmetrized_delta_nonzero': {f'{u},{v}': z for (u, v), z in sorted(sym.items())},
        'symmetrized_row_sums': {str(u): z for u, z in row_sums.items()},
        'symmetric_total': sum(sym.values()),
        'symmetric_e2_total': sum(geometry(u, v)[1] * z for (u, v), z in sym.items()),
        'nonzero_geometry_delta': {f'{d},{e}': z for (d, e), z in sorted(aggregate.items()) if z},
    }
    Path(__file__).with_name('n3_rootpair_flux_audit.json').write_text(
        json.dumps(out, sort_keys=True, indent=2) + '\n')
    print(json.dumps(out, sort_keys=True, indent=2))


if __name__ == '__main__':
    main()
