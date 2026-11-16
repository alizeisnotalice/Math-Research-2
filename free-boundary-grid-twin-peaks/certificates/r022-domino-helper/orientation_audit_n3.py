#!/usr/bin/env python3
"""Full 3x3 word audit of oriented-cover collapse in the K3 helper map."""
from collections import Counter, defaultdict
from itertools import permutations
import json

n = 3
vertices = tuple((x, y) for y in range(n) for x in range(n))
vset = set(vertices)
nbr = {v: {u for u in ((v[0] - 1, v[1]), (v[0] + 1, v[1]),
                        (v[0], v[1] - 1), (v[0], v[1] + 1)) if u in vset}
       for v in vertices}
edges = [(a, b) for a in vertices for b in nbr[a]]  # deliberately oriented


def peaks(word):
    seen, ans = set(), []
    for x in word:
        if not (nbr[x] & seen):
            ans.append(x)
        seen.add(x)
    return ans


def common(p):
    ans = set(vertices)
    for x in p:
        ans &= nbr[x]
    return ans


by_source_mark = defaultdict(list)
checked_targets = 0
k3_words = 0
oriented_records = 0
for word in permutations(vertices):
    pset = set(peaks(word))
    if len(pset) != 3:
        continue
    k3_words += 1
    if common(pset):
        continue
    pos = {x: i for i, x in enumerate(word)}
    for a, b in edges:
        ca, cb = nbr[a] - {b}, nbr[b] - {a}
        if not pset <= ({a} | ca | {b} | cb):
            continue
        ka, kb = len(pset & ca), len(pset & cb)
        if sorted((ka, kb)) != [1, 2]:
            continue
        # Normalize the heavy side to h; this deliberately forgets the
        # original orientation (a,b), which the audit measures.
        h, ell = (a, b) if ka == 2 else (b, a)
        heavy = sorted(pset & (nbr[h] - {ell}), key=pos.__getitem__)
        p, r = heavy
        q = next(iter(pset - {p, r}))
        assert pos[p] < pos[r] < pos[h]
        assert q not in nbr[h]
        target = list(word)
        target.remove(h)
        target.insert(target.index(p) + 1, h)
        assert set(peaks(target)) == {p, q}
        mark = (p, h, q, r, ell)
        by_source_mark[(word, mark, tuple(target))].append((a, b))
        oriented_records += 1
        checked_targets += 1

hist = Counter(len(v) for v in by_source_mark.values())
assert hist == Counter({2: len(by_source_mark)})
for orientations in by_source_mark.values():
    assert set(orientations) == {(orientations[0][0], orientations[0][1]),
                                 (orientations[0][1], orientations[0][0])}

print(json.dumps({
    'status': 'passed',
    'grid': '3x3',
    'all_K3_words': k3_words,
    'oriented_cover_helper_records': oriented_records,
    'normalized_source_mark_target_records': len(by_source_mark),
    'oriented_multiplicity_histogram': dict(sorted(hist.items())),
    'all_paired_orientations_have_identical_normalized_target': True,
    'checked_target_peaksets': checked_targets,
}, indent=2, sort_keys=True))
