#!/usr/bin/env python3
"""Independent-set inclusion-exclusion count of A1 for the 5x5 corner deletion.

This does not use connected-order subset DP.  After complementing labels,
one-peak labelings become labelings with exactly one local minimum.  For an
independent S, forcing every vertex of S to be a local minimum is the
height-two poset s < u for all edges su.  The identity

  1_{number of minima = 1} = sum_{empty != S subset minima} (-1)^(|S|-1)|S|

reduces A1 to independent-set inclusion-exclusion.
"""

from collections import Counter
from math import factorial
from pathlib import Path
import json
import time
import argparse


N = 5
D = {(0, 0), (1, 0), (0, 1)}
V = tuple((x, y) for y in range(N) for x in range(N) if (x, y) not in D)
IDX = {v: i for i, v in enumerate(V)}
ADJ = []
for x, y in V:
    z = 0
    for w in ((x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1)):
        j = IDX.get(w)
        if j is not None:
            z |= 1 << j
    ADJ.append(z)


def independent_masks():
    """All independent subsets of H, each exactly once, in lexicographic DFS."""
    m = len(V)

    def visit(i, chosen, forbidden):
        if i == m:
            yield chosen
            return
        yield from visit(i + 1, chosen, forbidden)
        if not (forbidden >> i & 1):
            yield from visit(i + 1, chosen | (1 << i), forbidden | ADJ[i] | (1 << i))

    yield from visit(0, 0, 0)


def linear_extensions_forced_minima(smask):
    """Exact extension count for selected independent minima (height two)."""
    selected = [i for i in range(len(V)) if smask >> i & 1]
    m = len(selected)
    if not m:
        return factorial(len(V))
    source_bit = {v: 1 << j for j, v in enumerate(selected)}
    need = [0] * len(V)
    for u in range(len(V)):
        a = ADJ[u] & smask
        while a:
            low = a & -a
            v = low.bit_length() - 1
            need[u] |= source_bit[v]
            a ^= low
    assert all(need[v] == 0 for v in selected)
    nonsources = [u for u in range(len(V)) if u not in source_bit]
    q = len(nonsources)
    # e[A] is # non-sources whose required selected minima are all in A.
    eligible = [sum((need[u] & ~a) == 0 for u in nonsources) for a in range(1 << m)]
    dp = [[0] * (q + 1) for _ in range(1 << m)]
    dp[0][0] = 1
    full = (1 << m) - 1
    for a in range(1 << m):
        for k, value in enumerate(dp[a]):
            if not value:
                continue
            unused = full ^ a
            while unused:
                low = unused & -unused
                unused ^= low
                dp[a | low][k] += value
            if k < q:
                choices = eligible[a] - k
                if choices:
                    dp[a][k + 1] += value * choices
    return dp[full][q]


def main(summary_path: Path):
    started = time.perf_counter()
    signed = 0
    hist = Counter()
    positive = negative = 0
    terms = 0
    for s in independent_masks():
        k = s.bit_count()
        hist[k] += 1
        if not k:
            continue
        le = linear_extensions_forced_minima(s)
        term = k * le
        terms += 1
        if k & 1:
            signed += term
            positive += term
        else:
            signed -= term
            negative += term
    # From global tail exact table: tail_T_ge_t[21] is T_22 because index 0=t=1.
    data = json.loads(summary_path.read_text())
    t22 = int(data['tail_T_ge_t'][21])
    result = {
        'method': 'independent-set inclusion-exclusion for exactly one local minimum after label complement',
        'n': 5,
        'H_vertices': len(V),
        'deleted_D': sorted(D),
        'independent_sets_including_empty': sum(hist.values()),
        'independent_size_histogram': {str(k): hist[k] for k in sorted(hist)},
        'nonempty_IE_terms': terms,
        'positive_sum': str(positive),
        'negative_sum': str(negative),
        'A1_H': str(signed),
        'global_source_sha256': __import__('hashlib').sha256(summary_path.read_bytes()).hexdigest(),
        'T_22_from_global_tail': str(t22),
        'eight_A1_H': str(8 * signed),
        'T_22_equals_eight_A1_H': t22 == 8 * signed,
        'elapsed_seconds': time.perf_counter() - started,
    }
    print(json.dumps(result, indent=2, sort_keys=True))
    assert t22 == 8 * signed


if __name__ == '__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--tail-summary',type=Path,required=True)
    main(ap.parse_args().tail_summary)
