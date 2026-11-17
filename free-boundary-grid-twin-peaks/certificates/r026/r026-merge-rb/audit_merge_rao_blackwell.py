#!/usr/bin/env python3
"""Exact n=3 and finite-trace audit for first-merge Rao--Blackwell weights.

The statistic m=q(q-1)/(2ab) is a conditional mean over an entire fixed
first-merge fibre.  This script deliberately distinguishes that identity from
the false pointwise assertion L=m for an individual tau=2 word.
"""
from __future__ import annotations

import hashlib
import os
import itertools
import json
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
CERT_ROOT = Path(os.environ.get('A022_CERTIFICATES_ROOT', HERE.parent.parent))
TRACE_DIR = CERT_ROOT / 'r025' / 'r025-rank-bias'
TRACE_RAW = TRACE_DIR / 'n3-p17-r1-l4-trace.jsonl'
TRACE_JSON = TRACE_DIR / 'n3-p17-r1-l4-trace.jsonl.trace.json'


def graph(n: int) -> list[set[int]]:
    ans = [set() for _ in range(n * n)]
    for v in range(n * n):
        x, y = v % n, v // n
        for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            xx, yy = x + dx, y + dy
            if 0 <= xx < n and 0 <= yy < n:
                ans[v].add(yy * n + xx)
    return ans


def components(adj: list[set[int]], vertices: set[int]) -> list[set[int]]:
    remaining, out = set(vertices), []
    while remaining:
        pending, comp = [remaining.pop()], set()
        while pending:
            v = pending.pop()
            comp.add(v)
            found = adj[v] & remaining
            remaining -= found
            pending.extend(found)
        out.append(comp)
    return out


def births(adj: list[set[int]], word: tuple[int, ...]) -> tuple[int, ...]:
    seen, out = set(), []
    for rank, v in enumerate(word, 1):
        if not (adj[v] & seen):
            out.append(rank)
        seen.add(v)
    return tuple(out)


def merge_data(adj: list[set[int]], word: tuple[int, ...]) -> tuple[tuple[object, ...], int, int, int, int]:
    """Return a complete fixed-fibre key and a,b,q,merge rank."""
    seen: set[int] = set()
    for rank, v in enumerate(word, 1):
        cc = components(adj, seen) if seen else []
        touched = [c for c in cc if adj[v] & c]
        if len(touched) >= 2:
            assert len(cc) == len(touched) == 2
            A, B = sorted(touched, key=lambda c: min(c))
            a, b = len(A), len(B)
            seqA = tuple(x for x in word[:rank - 1] if x in A)
            seqB = tuple(x for x in word[:rank - 1] if x in B)
            key = (tuple(sorted(word[:rank])), v, seqA, seqB, tuple(word[rank:]))
            return key, a, b, a + b, rank
        seen.add(v)
    raise AssertionError('missing first merge')


def merge_sizes(adj: list[set[int]], word: tuple[int, ...]) -> tuple[int, int, int]:
    _, a, b, q, _ = merge_data(adj, word)
    return a, b, q


def root_observables(n: int, word: tuple[int, ...], adj: list[set[int]]) -> tuple[int, int, int]:
    br = births(adj, word)
    assert len(br) == 2
    x, y = word[br[0] - 1], word[br[1] - 1]
    dx, dy = x % n - y % n, x // n - y // n
    d2 = dx * dx + dy * dy
    ax, ay = 2 * (x % n) - (n - 1), 2 * (x // n) - (n - 1)
    bx, by = 2 * (y % n) - (n - 1), 2 * (y // n) - (n - 1)
    f0, g0 = ax * ax + ay * ay + bx * bx + by * by, ax * bx + ay * by
    assert f0 - 2 * g0 == 4 * d2
    return d2, f0, g0


def observables(n: int, word: tuple[int, ...], adj: list[set[int]]) -> tuple[int, int, int, int, int]:
    assert births(adj, word) == (1, 2)
    _, y = word[:2]
    rank = next(k for k, v in enumerate(word, 1) if k >= 3 and v in adj[y])
    L = rank - 2
    d2, f0, g0 = root_observables(n, word, adj)
    return L, d2, f0, g0, rank


def add(dst: dict[str, Fraction], weight: Fraction, h: tuple[int, int, int, int]) -> None:
    dst['one'] += weight
    for name, value in zip(('d2', 'f0', 'g0'), h[1:4]):
        dst[name] += weight * value


def as_text(x: dict[str, Fraction]) -> dict[str, str]:
    return {k: str(v) for k, v in x.items()}


def n3_exact() -> dict[str, object]:
    n, N, adj = 3, 9, graph(3)
    full = defaultdict(Fraction)
    tau2_L = defaultdict(Fraction)
    tau2_m = defaultdict(Fraction)
    fibre_full: dict[tuple[object, ...], dict[str, Fraction]] = defaultdict(lambda: defaultdict(Fraction))
    fibre_tau2_L: dict[tuple[object, ...], dict[str, Fraction]] = defaultdict(lambda: defaultdict(Fraction))
    fibre_tau2_m: dict[tuple[object, ...], dict[str, Fraction]] = defaultdict(lambda: defaultdict(Fraction))
    first_merge_at_last = 0
    individual_inequality = 0
    for word in itertools.permutations(range(N)):
        br = births(adj, word)
        if len(br) != 2:
            continue
        d2, f0, g0 = root_observables(n, word, adj)
        fibre, a_all, b_all, q_all, merge_rank = merge_data(adj, word)
        for name, value in zip(('one', 'd2', 'f0', 'g0'), (1, d2, f0, g0)):
            full[name] += value
            fibre_full[fibre][name] += value
        if br[1] != 2:
            continue
        L, _, _, _, _ = observables(n, word, adj)
        h = (L, d2, f0, g0)
        a, b, q = a_all, b_all, q_all
        if merge_rank == N:
            first_merge_at_last += 1
        m = Fraction(q * (q - 1), 2 * a * b)
        if Fraction(L) != m:
            individual_inequality += 1
        add(tau2_L, Fraction(L), h)
        add(tau2_m, m, h)
        add(fibre_tau2_L[fibre], Fraction(L), h)
        add(fibre_tau2_m[fibre], m, h)
    assert full == tau2_L == tau2_m
    assert set(fibre_full) == set(fibre_tau2_L) == set(fibre_tau2_m)
    assert all(fibre_full[k] == fibre_tau2_L[k] == fibre_tau2_m[k] for k in fibre_full)
    assert full['one'] == 143808 and full['d2'] == 633600
    return {
        'method': 'all 9! words, exact Fraction arithmetic, independent first-merge scan',
        'full_K2_totals': as_text(full),
        'tau2_L_totals': as_text(tau2_L),
        'tau2_m_totals': as_text(tau2_m),
        'every_target_total_recovers_exactly': True,
        'complete_fixed_fibres': len(fibre_full),
        'every_complete_fibre_recovers_all_H_totals': True,
        'tau2_words_with_L_not_equal_m': individual_inequality,
        'tau2_words_first_merge_at_terminal_last_vertex': first_merge_at_last,
    }


def trace_diagnostic() -> dict[str, object]:
    raw_lines = [json.loads(x) for x in TRACE_RAW.read_text().splitlines() if x.strip()]
    trace = json.loads(TRACE_JSON.read_text())
    meta, rep = raw_lines[0], next(x for x in raw_lines if x['kind'] == 'replicate')
    assert meta['implementation'] == 'r025-tau2-rankbiased-prefix-mh'
    terminal = next(x for x in trace['rows'] if x.get('terminal'))
    n, adj = trace['n'], graph(trace['n'])
    Lweighted = defaultdict(Fraction)
    mweighted = defaultdict(Fraction)
    last_merge = 0
    pointwise_unequal = 0
    records = []
    for prefix in terminal['prefix']:
        missing = next(v for v in range(n * n) if v not in prefix)
        word = tuple(prefix + [missing])
        L, d2, f0, g0, _ = observables(n, word, adj)
        a, b, q = merge_sizes(adj, word)
        m = Fraction(q * (q - 1), 2 * a * b)
        pointwise_unequal += Fraction(L) != m
        add(Lweighted, Fraction(L), (L, d2, f0, g0))
        add(mweighted, m, (L, d2, f0, g0))
        # This independently recognizes whether the terminal append was the
        # first merge; the online design must therefore scan the omitted point.
        seen = set(word[:-1])
        cc = components(adj, seen)
        if len([c for c in cc if adj[word[-1]] & c]) >= 2:
            last_merge += 1
        records.append({'a': a, 'b': b, 'q': q, 'm': str(m), 'L': L, 'd2': d2})
    P = trace['P']
    assert int(rep['tau2_lift_L_groups'][0]) == Lweighted['one']
    assert int(rep['tau2_lift_L_distance2_groups'][0]) == Lweighted['d2']
    assert int(rep['tau2_lift_L_F0_groups'][0]) == Lweighted['f0']
    assert int(rep['tau2_lift_L_G0_groups'][0]) == Lweighted['g0']
    assert Lweighted['f0'] - 2 * Lweighted['g0'] == 4 * Lweighted['d2']
    assert mweighted['f0'] - 2 * mweighted['g0'] == 4 * mweighted['d2']
    return {
        'input_raw_sha256': hashlib.sha256(TRACE_RAW.read_bytes()).hexdigest(),
        'input_trace_sha256': hashlib.sha256(TRACE_JSON.read_bytes()).hexdigest(),
        'particles': P,
        'particles_with_L_not_equal_m': pointwise_unequal,
        'L_weighted_terminal_totals': as_text(Lweighted),
        'm_weighted_terminal_totals': as_text(mweighted),
        'terminal_last_vertex_is_first_merge_for_particles': last_merge,
        'per_particle_first_merge': records,
        'scope': 'Finite trace comparison only. Different L and m totals on this incomplete particle sample do not contradict the fibre conditional-expectation identity.',
    }


def main() -> None:
    output = {
        'scope': 'Exact target-law n3 identity plus finite r025 terminal-trace diagnostic; no GPU execution',
        'n3_exact': n3_exact(),
        'trace_diagnostic': trace_diagnostic(),
        'theory': 'For a complete fixed first-merge fibre and a root-symmetric interleaving-invariant H, E_tau2[L H | fibre] = m H. This is not L=m per word and does not itself establish finite-SMC variance reduction.',
    }
    HERE.joinpath('merge-rb-audit.json').write_text(json.dumps(output, indent=2) + '\n')
    print(json.dumps(output, indent=2))


if __name__ == '__main__':
    main()
