#!/usr/bin/env python3
"""Exact support-tree classification of pinned R63 global maximum blocks.

Reads saved input/max-mask/component metadata only; no halo or cut evaluation.
The saved eligible-component partition is inherited, with its masses and span
ratios checked afresh. Global maximum masks are inherited, not re-enumerated.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from functools import cmp_to_key
import hashlib
import json
from pathlib import Path
import platform
import time


def radical_sign(a, b, L):
    """Sign of a*sqrt(L)-b, using rational comparisons only."""
    if not a:
        return (b < 0) - (b > 0)
    if a < 0:
        return -radical_sign(-a, -b, L)
    if b <= 0:
        return 1
    v = L*a*a-b*b
    return (v > 0)-(v < 0)


def mst(points, ids, L):
    edges = []
    for k, i in enumerate(ids):
        for j in ids[k+1:]:
            dx, dy = (abs(points[i][t]-points[j][t]) for t in (0, 1))
            edges.append((dx+dy, dx*dy, i, j))
    def compare(e, f):
        z = radical_sign(e[0]-f[0], e[1]-f[1], L)
        return z or ((e[2:] > f[2:])-(e[2:] < f[2:]))
    edges.sort(key=cmp_to_key(compare))
    parent = {i:i for i in ids}
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    used, P, Q = [], F(0), F(0)
    for p, q, i, j in edges:
        a, b = find(i), find(j)
        if a != b:
            parent[a] = b
            used.append([i, j]); P += p; Q += q
    assert len(used) == len(ids)-1
    target = F(5, 2)*L+Q
    margin = target*target-L*P*P
    return {'edges':used, 'P':str(P), 'Q':str(Q),
            'squared_test_margin':str(margin), 'paid':margin >= 0}


def mask_ids(mask, n):
    return [i for i in range(n) if mask >> i & 1]


def geometry(mask, points, masses):
    ids = mask_ids(mask, len(points))
    mass = sum((masses[i] for i in ids), F(0))
    widths = [max(points[i][t] for i in ids)-min(points[i][t] for i in ids) for t in (0,1)]
    return ids, mass, widths


def counts(rows):
    return {'blocks':len(rows), 'priority_classes':dict(sorted(Counter(r['priority_class'] for r in rows).items())),
            'overlapping_flags':{k:sum(r[k] for r in rows) for k in
              ['component_terminal','whole_input_terminal','thin','rank_at_most_3','tree_paid','proper']}}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--archive-root', type=Path, required=True, help='Topic directory containing archived computation/results/round063.')
    ap.add_argument('--spec', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args(); started = time.perf_counter()
    spec = json.loads(args.spec.read_text())
    rows, case_rows = [], []
    for src in spec['sources']:
        path = args.archive_root / src['path']
        blob = path.read_bytes()
        assert len(blob) == src['bytes'] and hashlib.sha256(blob).hexdigest() == src['sha256']
        data = json.loads(blob)
        assert data['status'] == 'PASS'
        for case in data['cases']:
            inp = case['input']; n = len(inp['points_quarters']); full = (1 << n)-1
            assert inp['coordinate_unit'] == '1/4' and inp['mass_unit'] == '1/16'
            pts = [tuple(F(x,4) for x in xy) for xy in inp['points_quarters']]
            masses = [F(m,16) for m in inp['mass_sixteenths']]
            assert all(m > 0 for m in masses) and len(set(pts)) == n
            _, total, fw = geometry(full,pts,masses)
            full_terminal = max(fw)**2 <= F(35,33)*total
            comps = case['eligible_components']; partition = 0
            for comp in comps:
                cm = comp['mask']; assert not partition & cm; partition |= cm
                ci, mass, widths = geometry(cm,pts,masses)
                assert len(ci) == comp['rank']
                assert max(widths)**2/mass == F(comp['span_squared_over_mass'])
                assert (max(widths)**2 <= F(35,33)*mass) == comp['terminal_35_33']
            assert partition == full
            de = case['deletions']; L = F(de['maximum_strictly_eligible_mass_sixteenths'],16)
            masks = de['all_maximum_strictly_eligible_masks']
            assert len(set(masks)) == len(masks)
            current = []
            for mask in masks:
                ids, mass, widths = geometry(mask,pts,masses)
                assert mass == L and max(widths)**2 < L
                # Verify bbox closure of every inherited maximum mask.
                low = [min(pts[i][t] for i in ids) for t in (0,1)]
                high = [max(pts[i][t] for i in ids) for t in (0,1)]
                closure = sum(1<<i for i,p in enumerate(pts) if all(low[t] <= p[t] <= high[t] for t in (0,1)))
                assert closure == mask
                containing = [c for c in comps if mask & c['mask'] == mask]
                assert len(containing) == 1
                ct = containing[0]['terminal_35_33']
                thin = min(widths)**2 <= F(9,16)*L
                small = len(ids) <= 3
                tree = mst(pts,ids,L)
                category = ('component_terminal' if ct else 'thin' if thin else
                            'rank_at_most_3' if small else 'support_tree' if tree['paid'] else 'remaining')
                row = {'input_id':inp['id'], 'source_count':n, 'mask':mask, 'rank':len(ids),
                       'L':str(L), 'widths':list(map(str,widths)),
                       'width_squared_over_L':[str(v*v/L) for v in widths],
                       'proper':mask != full, 'eligible_component_mask':containing[0]['mask'],
                       'eligible_component_count':len(comps), 'component_terminal':ct,
                       'whole_input_terminal':full_terminal, 'thin':thin,
                       'rank_at_most_3':small, 'tree_paid':tree['paid'],
                       'priority_class':category, 'mst':tree}
                current.append(row); rows.append(row)
            assert bool(full in masks) == de['maximum_strictly_eligible_full_input']
            case_rows.append({'input_id':inp['id'], 'maximum_ties':len(current),
                              'connected':len(comps)==1, 'whole_input_terminal':full_terminal,
                              'some_maximum_old_branch_paid':any(r['component_terminal'] or r['thin'] or r['rank_at_most_3'] for r in current),
                              'some_maximum_tree_paid':any(r['tree_paid'] for r in current),
                              'some_maximum_classified_paid':any(r['priority_class'] != 'remaining' for r in current),
                              'all_maximum_unclassified':all(r['priority_class'] == 'remaining' for r in current)})
    assert len(rows) == 211 and len(case_rows) == 72 and sum(r['proper'] for r in rows) == 199
    remaining = [r for r in rows if r['priority_class'] == 'remaining']
    summary = {'all_global_maximum_ties':counts(rows), 'proper_global_maximum_ties':counts([r for r in rows if r['proper']]),
               'case_count':len(case_rows), 'case_flags':{k:sum(r[k] for r in case_rows) for k in
                  ['connected','whole_input_terminal','some_maximum_old_branch_paid','some_maximum_tree_paid','some_maximum_classified_paid','all_maximum_unclassified']},
               'remaining_rows':len(remaining), 'remaining_input_ids':sorted(set(r['input_id'] for r in remaining))}
    result = {'schema':'r068-archived-global-max-tree-v1', 'status':'PASS',
              'scope':'211 inherited global maximum ties, not all component maxima. Exact Kruskal edge ordering and exact tree payment. No halo/cut replay and no new selection theorem.',
              'inherited_metadata':'The R63 maximum masks and eligible-component partitions are hash-bound; this classifier verifies their geometry but does not independently regenerate those catalogues.',
              'python':platform.python_version(), 'platform':platform.platform(),
              'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'spec_sha256':hashlib.sha256(args.spec.read_bytes()).hexdigest(),
              'summary':summary, 'cases':case_rows, 'blocks':rows,
              'wall_seconds':time.perf_counter()-started}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'PASS','summary':summary,'wall_seconds':result['wall_seconds']}))


if __name__ == '__main__':
    main()
