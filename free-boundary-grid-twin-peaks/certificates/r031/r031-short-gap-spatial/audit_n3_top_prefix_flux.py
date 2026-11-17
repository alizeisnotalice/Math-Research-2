#!/usr/bin/env python3
"""Independent n=3 check of the top-four/top-five completion formula."""

from collections import defaultdict
import json
from pathlib import Path


n = 3
N = n * n
V = tuple(range(N))
ADJ = [0] * N
for v in V:
    x, y = v % n, v // n
    for xx, yy in ((x + 1, y), (x - 1, y), (x, y + 1), (x, y - 1)):
        if 0 <= xx < n and 0 <= yy < n:
            ADJ[v] |= 1 << (yy * n + xx)
FULL = (1 << N) - 1


def completion_count():
    f = [0] * (FULL + 1)
    f[FULL] = 1
    for S in range(FULL - 1, 0, -1):
        outside = FULL ^ S
        while outside:
            bit = outside & -outside
            outside ^= bit
            x = bit.bit_length() - 1
            if ADJ[x] & S:
                f[S] += f[S | bit]
    return f


def main():
    f = completion_count()
    h2, h3 = defaultdict(int), defaultdict(int)
    for v in V:
        bv = 1 << v
        for p in V:
            bp = 1 << p
            if p == v or ADJ[v] & bp:
                continue
            for q in V:
                bq = 1 << q
                if q in (v, p) or not (ADJ[v] & bq) or ADJ[p] & bq:
                    continue
                for x in V:
                    bx = 1 << x
                    if x in (v, p, q) or not (ADJ[p] & bx):
                        continue
                    h2[v, p] += f[bv | bp | bq | bx]
            for q in V:
                bq = 1 << q
                if q in (v, p) or not (ADJ[v] & bq) or ADJ[p] & bq:
                    continue
                for r in V:
                    br = 1 << r
                    if r in (v, p, q) or ADJ[p] & br or not (ADJ[r] & (bv | bq)):
                        continue
                    for x in V:
                        bx = 1 << x
                        if x in (v, p, q, r) or not (ADJ[p] & bx):
                            continue
                        h3[v, p] += f[bv | bp | bq | br | bx]
    prior = json.loads(Path(__file__).with_name('n3_rootpair_flux_audit.json').read_text())
    ref2 = {tuple(map(int, s.split(','))): int(z) for s, z in prior['directed_h2'].items()}
    ref3 = {tuple(map(int, s.split(','))): int(z) for s, z in prior['directed_h3'].items()}
    assert dict(h2) == ref2
    assert dict(h3) == ref3
    out = {
        'status': 'passed',
        'h2_total': str(sum(h2.values())),
        'h3_total': str(sum(h3.values())),
        'formula_checked': 'top-four/top-five prefix sums against completion f_G(S)',
    }
    Path(__file__).with_name('n3_top_prefix_flux_audit.json').write_text(
        json.dumps(out, sort_keys=True, indent=2) + '\n')
    print(json.dumps(out, sort_keys=True, indent=2))


if __name__ == '__main__':
    main()
