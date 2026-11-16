#!/usr/bin/env python3
"""Exact reflection-pair orientation orbits for central vertical reflection.

For fixed unordered labels on each reflected vertex pair, enumerate every
orientation.  Count only words/labelings with exactly two peaks and classify
their peaks as same-side or split.  n=2 is exhaustive over all pair orbits;
n=4 is a deterministic random-orbit search for a strict counterexample.
"""
from itertools import combinations
import json
import random


def reflection_pairs(n):
    return [((x, y), (n - 1 - x, y))
            for y in range(n) for x in range(n // 2)]


def orbit_counts(n, label_pairs):
    pairs = reflection_pairs(n)
    assert len(pairs) == len(label_pairs) == n * n // 2
    side = {(x, y): x < n // 2 for y in range(n) for x in range(n)}
    neighbours = {(x, y): [(xx, yy) for xx, yy in
                            ((x - 1, y), (x + 1, y),
                             (x, y - 1), (x, y + 1))
                            if 0 <= xx < n and 0 <= yy < n]
                  for y in range(n) for x in range(n)}
    same = split = 0
    peak_hist = {}
    for spins in range(1 << len(pairs)):
        lab = {}
        for i, (left, right) in enumerate(pairs):
            a, b = label_pairs[i]
            if (spins >> i) & 1:
                a, b = b, a
            lab[left], lab[right] = a, b
        peaks = [v for v in lab if all(lab[v] > lab[u] for u in neighbours[v])]
        peak_hist[len(peaks)] = peak_hist.get(len(peaks), 0) + 1
        if len(peaks) == 2:
            if side[peaks[0]] == side[peaks[1]]:
                same += 1
            else:
                split += 1
    return {'same': same, 'split': split, 'peak_histogram': peak_hist}


def two_peak_orientations(n, label_pairs):
    """Compact certificate: spin bit i means swap the two labels in pair i."""
    pairs = reflection_pairs(n)
    side = {(x, y): x < n // 2 for y in range(n) for x in range(n)}
    neighbours = {(x, y): [(xx, yy) for xx, yy in
                            ((x - 1, y), (x + 1, y),
                             (x, y - 1), (x, y + 1))
                            if 0 <= xx < n and 0 <= yy < n]
                  for y in range(n) for x in range(n)}
    out = []
    for spins in range(1 << len(pairs)):
        lab = {}
        for i, (left, right) in enumerate(pairs):
            a, b = label_pairs[i]
            if (spins >> i) & 1: a, b = b, a
            lab[left], lab[right] = a, b
        peaks = [v for v in lab if all(lab[v] > lab[u] for u in neighbours[v])]
        if len(peaks) == 2:
            out.append({'spins_low_to_high_pair_order': format(spins, f'0{len(pairs)}b'),
                        'peaks': peaks,
                        'class': 'same' if side[peaks[0]] == side[peaks[1]] else 'split'})
    return out


def exhaustive_n2():
    # The two geometric reflection pairs are distinguished by their row.
    out = []
    for a, b in combinations(range(1, 5), 2):
        rem = sorted(set(range(1, 5)) - {a, b})
        pairs = [(a, b), tuple(rem)]
        c = orbit_counts(2, pairs)
        out.append({'label_pairs': pairs, **c})
    return out


def search_n4(limit=20000, seed=20261116):
    rng = random.Random(seed)
    labels = list(range(1, 17))
    best = None
    for trial in range(limit):
        rng.shuffle(labels)
        pairs = [tuple(sorted(labels[2*i:2*i+2])) for i in range(8)]
        c = orbit_counts(4, pairs)
        record = {'trial': trial, 'label_pairs': pairs, **c}
        if best is None or c['same'] - c['split'] > best['same'] - best['split']:
            best = record
        if c['same'] > c['split']:
            record['two_peak_orientations'] = two_peak_orientations(4, pairs)
            return {'found': True, 'searched': trial + 1, 'witness': record,
                    'best': best}
    return {'found': False, 'searched': limit, 'best': best}


n2 = exhaustive_n2()
answer = {
    'status': 'passed',
    'n2_orbits': n2,
    'n2_all_split_ge_same': all(x['split'] >= x['same'] for x in n2),
    'n4_search': search_n4(),
}
print(json.dumps(answer, indent=2, sort_keys=True))
