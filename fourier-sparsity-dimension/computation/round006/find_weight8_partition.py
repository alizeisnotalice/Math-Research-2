#!/usr/bin/env python3
"""Find a weight-8 exact-cover partition inside each dim-15 sign-code V.

Rows with a common 3-space U receive its eight cosets directly: dimension 15
equals the affine-on-U-cosets parameter-space bound.  For the remaining rows,
enumerate all 2^15 codewords in C plus the eight stored generators and solve
an exact-cover problem using only weight-8 codewords.
"""
import argparse
import hashlib
import json
from pathlib import Path

POINTS, FULL = 64, (1 << 64) - 1

def rank(vectors):
    piv = {}
    for x in vectors:
        while x:
            b = x.bit_length() - 1
            if b in piv: x ^= piv[b]
            else: piv[b] = x; break
    return len(piv)

def c_basis():
    return [FULL] + [sum(((x >> j) & 1) << x for x in range(POINTS)) for j in range(6)]

def coset_masks(basis):
    U = {0}
    for b in basis: U |= {u ^ b for u in list(U)}
    answer, seen = [], set()
    for x in range(POINTS):
        m = sum(1 << (x ^ u) for u in U)
        if m not in seen: seen.add(m); answer.append(m)
    if len(answer) != 8 or any(m.bit_count() != 8 for m in answer): raise AssertionError(basis)
    return sorted(answer)

def enumerate_weight8(generators):
    if len(generators) != 15 or rank(generators) != 15: raise AssertionError('V basis')
    candidates = {}
    for z in range(1 << 15):
        m = 0
        for j, g in enumerate(generators):
            if z >> j & 1: m ^= g
        if m.bit_count() == 8:
            candidates[m] = z
    return candidates

def first_exact_cover(candidates):
    by_point = [[] for _ in range(POINTS)]
    for mask in candidates:
        for x in range(POINTS):
            if mask >> x & 1: by_point[x].append(mask)
    def dfs(remaining, chosen):
        if not remaining: return chosen
        choices = None
        points = remaining
        while points:
            low = points & -points; x = low.bit_length() - 1; points ^= low
            legal = [m for m in by_point[x] if not (m & ~remaining)]
            if not legal: return None
            if choices is None or len(legal) < len(choices):
                choices = legal
                if len(choices) == 1: break
        for mask in choices:
            got = dfs(remaining ^ mask, chosen + [mask])
            if got is not None: return got
        return None
    return dfs(FULL, [])

parser = argparse.ArgumentParser()
parser.add_argument('audit', type=Path)
parser.add_argument('--out', type=Path, required=True)
args = parser.parse_args()
raw = args.audit.read_bytes(); audit = json.loads(raw)
if audit['certificate_count'] != 20 or audit['subspaces_tested'] != 1395: raise AssertionError('audit scope')
rows = []
for row in audit['rows']:
    generators = c_basis() + row['relative_generators']
    if row['valid_U_count']:
        blocks = coset_masks(row['first_valid_U_rref_basis'])
        rows.append({'source_file': row['source_file'], 'row_index': row['row_index'], 'seed_index': row['seed_index'],
                     'method': 'common_U_cosets_by_dim15_equality', 'weight8_candidate_count': None,
                     'partition_found': True, 'partition_masks_hex': [hex(x) for x in blocks],
                     'partition_codewords_hex': None, 'common_U_rref_basis': row['first_valid_U_rref_basis']})
    else:
        candidates = enumerate_weight8(generators)
        cover = first_exact_cover(candidates)
        rows.append({'source_file': row['source_file'], 'row_index': row['row_index'], 'seed_index': row['seed_index'],
                     'method': 'enumerate_2^15_weight8_then_exact_cover', 'weight8_candidate_count': len(candidates),
                     'partition_found': cover is not None,
                     'partition_masks_hex': [hex(x) for x in cover] if cover else None,
                     'partition_codewords_hex': [hex(candidates[x]) for x in cover] if cover else None,
                     'common_U_rref_basis': None})
if any(r['partition_found'] and len(r['partition_masks_hex']) != 8 for r in rows): raise AssertionError('block count')
args.out.write_text(json.dumps({'status': 'passed', 'audit_sha256': hashlib.sha256(raw).hexdigest(),
                                'method': 'common_U_or_full_2^15_weight8_exact_cover', 'rows': rows}, indent=2)+'\n')
print(json.dumps({'status':'passed','rows':len(rows),'partitions_found':sum(r['partition_found'] for r in rows),
                  'no_common_U_weight8_counts':[r['weight8_candidate_count'] for r in rows if r['method'].startswith('enumerate')]}))
