#!/usr/bin/env python3
"""Exact all-slot weighted accounting for stable helper moves on 3x3."""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import permutations
import json

n = 3
N = n*n
adj = []
for v in range(N):
    x, y = v % n, v // n
    adj.append({yy*n + xx for xx, yy in
                ((x-1,y), (x+1,y), (x,y-1), (x,y+1))
                if 0 <= xx < n and 0 <= yy < n})


def peak_list(word):
    seen, out = set(), []
    for v in word:
        if not (adj[v] & seen):
            out.append(v)
        seen.add(v)
    return out


def d(a, b):
    return abs(a % n-b % n) + abs(a//n-b//n)


load_num = Counter()
load_den = {}
quadruples = 0
source_words_with_quad = set()
weighted_total = Fraction(0)
k2d3_words = 0

for word_tuple in permutations(range(N)):
    word = list(word_tuple)
    roots = peak_list(word)
    if len(roots) == 2 and d(*roots) == 3:
        k2d3_words += 1
    if len(roots) != 3:
        continue
    pos = {v: i for i, v in enumerate(word)}
    rootset = set(roots)
    for p in roots:
        for q in roots:
            if p == q or d(p, q) != 3:
                continue
            r = next(iter(rootset - {p, q}))
            P, R = pos[p], pos[r]
            if P >= R:
                continue
            for z in range(N):
                Z = pos[z]
                if z in rootset or Z <= R:
                    continue
                if not (p in adj[z] and r in adj[z] and q not in adj[z]):
                    continue
                quadruples += 1
                source_words_with_quad.add(word_tuple)
                # There are g target slots after p and before r.  The reverse
                # old-z slot count is h=N-R-1 for every resulting target.
                g, h = R-P, N-R-1
                assert g > 0 and h > 0
                weighted_total += Fraction(g, h)
                for insertion in range(P + 1, R + 1):
                    new = word.copy()
                    new.pop(Z)          # Z>R, so r's index is unchanged.
                    new.insert(insertion, z)
                    assert set(peak_list(new)) == {p, q}
                    key = (tuple(new), p, z, q, r)
                    # In this target r lies one position later than in source.
                    target_h = N - new.index(r)
                    assert target_h == h
                    load_num[key] += 1
                    load_den[key] = h

assert all(v <= load_den[key] for key, v in load_num.items())
target_load = sum((Fraction(v, load_den[key]) for key, v in load_num.items()), Fraction(0))
assert target_load == weighted_total
assert all(d(key[1], key[3]) == 3 for key in load_num)

def frac(x):
    return {'numerator': x.numerator, 'denominator': x.denominator,
            'decimal': float(x)}

print(json.dumps({
    'status': 'passed', 'n': n,
    'K2_distance3_words': k2d3_words,
    'legal_ordered_helper_quadruples': quadruples,
    'K3_words_with_at_least_one_legal_quadruple': len(source_words_with_quad),
    'weighted_source_gain': frac(weighted_total),
    'weighted_gain_per_legal_quadruple': frac(weighted_total / quadruples),
    'weighted_gain_per_source_word_with_quadruple': frac(weighted_total / len(source_words_with_quad)),
    'marked_target_count': len(load_num),
    'max_preimage_count_per_marked_target': max(load_num.values(), default=0),
    'max_marked_target_load': frac(max((Fraction(v, load_den[k]) for k, v in load_num.items()), default=Fraction(0))),
    'preimage_count_histogram': {str(k): sum(v == k for v in load_num.values())
                                 for k in sorted(set(load_num.values()))},
}, indent=2, sort_keys=True))
