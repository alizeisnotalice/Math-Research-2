#!/usr/bin/env python3
"""Exact full-permutation audit of the added-nonedge response on P3 square."""
from collections import Counter, defaultdict
from itertools import combinations, permutations
import json

n = 3
vertices = tuple((x, y) for y in range(n) for x in range(n))
index = {v: i for i, v in enumerate(vertices)}
adj = [0] * len(vertices)
for i, (x, y) in enumerate(vertices):
    for j, (u, v) in enumerate(vertices):
        if abs(x - u) + abs(y - v) == 1:
            adj[i] |= 1 << j

nonedges = [(i, j) for i, j in combinations(range(len(vertices)), 2)
            if not (adj[i] >> j) & 1]
pair_count = Counter()
rooted_high = Counter()
base_a1 = 0
base_root_a1 = Counter()
direct_a1_plus = Counter()
direct_root_a1_plus = Counter()
response_condition_count = Counter()
off_endpoint_witness = {}


def peaks(word, changed=None):
    """Direct high-to-low peak evaluation; changed is an optional added edge."""
    seen, ans = 0, []
    x0 = y0 = None
    if changed is not None:
        x0, y0 = changed
    for v in word:
        nv = adj[v]
        if changed is not None:
            if v == x0:
                nv |= 1 << y0
            elif v == y0:
                nv |= 1 << x0
        if not (seen & nv):
            ans.append(v)
        seen |= 1 << v
    return ans


for word in permutations(range(len(vertices))):
    p = peaks(word)
    pos = {v: i for i, v in enumerate(word)}
    if len(p) == 1:
        base_a1 += 1
        base_root_a1[p[0]] += 1
    elif len(p) == 2:
        x, y = sorted(p)
        pair_count[(x, y)] += 1
        # Exactly one of the two roots is higher in the word.
        high = min(p, key=pos.__getitem__)
        low = y if high == x else x
        rooted_high[(high, low)] += 1
    # Independent direct check on every graph G+xy; it does not invoke the
    # claimed response identity.
    for e in nonedges:
        p_plus = peaks(word, e)
        if len(p_plus) == 1:
            direct_a1_plus[e] += 1
            direct_root_a1_plus[(e, p_plus[0])] += 1
        # This is the exact pre-addition event which makes the added edge
        # reduce the peak count.  The other endpoint need not be a peak.
        low = e[0] if pos[e[0]] > pos[e[1]] else e[1]
        if len(p) == 2 and low in p:
            response_condition_count[e] += 1
            if set(p) != set(e) and e not in off_endpoint_witness:
                off_endpoint_witness[e] = {
                    'word_high_to_low': [list(vertices[v]) for v in word],
                    'G_peakset': [list(vertices[v]) for v in p],
                    'lower_endpoint': list(vertices[low]),
                    'added_edge': [list(vertices[e[0]]), list(vertices[e[1]])],
                    'G_plus_edge_peakset': [list(vertices[v]) for v in peaks(word, e)],
                }

assert sum(pair_count.values()) == 143808
assert base_a1 == 29568
assert all(response_condition_count[e] == direct_a1_plus[e] - base_a1 for e in nonedges)
assert all(rooted_high[(x, y)] + rooted_high[(y, x)] == pair_count[(min(x, y), max(x, y))]
           for x, y in nonedges)
rooted_mismatches = []
for x, y in nonedges:
    for root, other in ((x, y), (y, x)):
        if direct_root_a1_plus[(x, y), root] - base_root_a1[root] != rooted_high[(root, other)]:
            rooted_mismatches.append((x, y, root, other,
                                      direct_root_a1_plus[(x, y), root] - base_root_a1[root],
                                      rooted_high[(root, other)]))
if rooted_mismatches:
    raise AssertionError({'first_rooted_mismatch': rooted_mismatches[0],
                          'coordinates': [vertices[z] for z in rooted_mismatches[0][:4]]})
assert all((len(vertices) - 1) * rooted_high[(x, y)] >= base_root_a1[x]
           for x, y in nonedges)

def dist(e):
    x, y = (vertices[z] for z in e)
    return abs(x[0] - y[0]) + abs(x[1] - y[1])

rows = []
for e in nonedges:
    x, y = e
    rows.append({
        'pair': [list(vertices[x]), list(vertices[y])],
        'distance': dist(e),
        'response': response_condition_count[e],
        'two_endpoint_peak_subcount': pair_count[e],
        'higher_root_counts': {
            f'{vertices[x]}>{vertices[y]}': rooted_high[(x, y)],
            f'{vertices[y]}>{vertices[x]}': rooted_high[(y, x)],
        },
    })
by_d = defaultdict(list)
for row in rows:
    by_d[row['distance']].append(row['response'])
summary = {str(d): {'pairs': len(v), 'minimum': min(v), 'maximum': max(v),
                    'total': sum(v)} for d, v in sorted(by_d.items())}

# Explicit witnesses that response is neither pointwise increasing nor
# pointwise decreasing with Manhattan length, even on P3 square.
inc_fail = next(((r, s) for r in rows for s in rows
                 if r['distance'] < s['distance'] and r['response'] > s['response']), None)
dec_fail = next(((r, s) for r in rows for s in rows
                 if r['distance'] < s['distance'] and r['response'] < s['response']), None)

print(json.dumps({
    'status': 'passed', 'grid': '3x3', 'A1_grid': base_a1,
    'rooted_A1_grid': {str(vertices[x]): base_root_a1[x] for x in range(len(vertices))},
    'nonedges': len(nonedges), 'response_distance_summary': summary,
    'pointwise_increasing_distance_failure': None if inc_fail is None else
        {'near': inc_fail[0], 'far': inc_fail[1]},
    'pointwise_decreasing_distance_failure': None if dec_fail is None else
        {'near': dec_fail[0], 'far': dec_fail[1]},
    'off_endpoint_witness_for_first_nonedge': off_endpoint_witness[nonedges[0]],
    'rooted_response_injection_minimum_slack': min(
        (len(vertices) - 1) * rooted_high[(x, y)] - base_root_a1[x]
        for x, y in nonedges),
    'pairs': rows,
}, indent=2, sort_keys=True))
