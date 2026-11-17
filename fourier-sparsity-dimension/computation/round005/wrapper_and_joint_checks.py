#!/usr/bin/env python3
"""Exact one-heavy-direction embeddings and failure of pairwise sign feasibility."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import platform
import sys


def walsh(v):
    v = v[:]
    step = 1
    while step < len(v):
        for start in range(0, len(v), 2 * step):
            for j in range(start, start + step):
                a, b = v[j], v[j + step]
                v[j], v[j + step] = a + b, a - b
        step *= 2
    return v


def rank(vectors):
    b = {}
    for x in vectors:
        while x:
            p = x.bit_length() - 1
            if p in b:
                x ^= b[p]
            else:
                b[p] = x
                break
    return len(b)


def pairs(support):
    return Counter(x ^ y for j, x in enumerate(support) for y in support[:j])


def check_wrapper(row):
    n, support, coefficients = row['n'], row['support'], row['integercoeffs']
    N = 1 << n
    char_bit, root_bit = N, 2 * N
    spectrum = {0: N, root_bit: N}
    for a, c in zip(support, coefficients):
        spectrum[a | char_bit] = c
        spectrum[a | char_bit | root_bit] = -c
    assert len(spectrum) == 2 * len(support) + 2
    frequency_list = sorted(spectrum)
    counts = pairs(frequency_list)
    assert counts[root_bit] == len(support) + 1
    assert all(counts[a] == 2 for a in frequency_list if a not in (0, root_bit))
    affine_rank_seed = rank([a ^ support[0] for a in support])
    assert rank(frequency_list) == affine_rank_seed + 2
    v = [spectrum.get(a, 0) for a in range(4 * N)]
    assert all(abs(x) == 2 * N for x in walsh(v))
    return len(support) > 1


def check_existing_normal_forms(row):
    n, support, coefficients = row['n'], row['support'], row['integercoeffs']
    if len(support) <= 4:
        return 0
    N = 1 << n
    coeff = dict(zip(support, coefficients))
    counts = pairs(support)
    checked = 0
    for anchor in support:
        heavy = [b for b in support if b != anchor and counts[anchor ^ b] >= 3]
        assert heavy
        if len(heavy) != 1:
            continue
        delta = anchor ^ heavy[0]
        assert {a ^ delta for a in support} == set(support)
        assert counts[delta] == len(support) // 2
        shifted = {a ^ anchor: c for a, c in coeff.items()}
        # Frequency translation is multiplication by a character. Physical
        # translation and overall sign make both anchor coefficients positive.
        u = 0 if shifted[0] * shifted[delta] > 0 else delta & -delta
        shifted = {a: c * (-1 if (a & u).bit_count() & 1 else 1)
                   for a, c in shifted.items()}
        sign = 1 if shifted[0] > 0 else -1
        shifted = {a: sign * c for a, c in shifted.items()}
        assert shifted[0] == shifted[delta] == N // 2
        assert all(shifted[a] == -shifted[a ^ delta] for a in shifted if a not in (0, delta))
        # Canonical quotient representatives use the pivot coordinate of delta.
        pivot = delta & -delta
        T = [a for a in shifted if a not in (0, delta) and not (a & pivot)]
        Tset = set(T)
        assert 0 not in Tset and not any((a ^ b) in Tset for a in T for b in T)
        checked += 1
    return checked


def joint_counterexample():
    rows = []
    for z in range(8):
        raw = []
        for x in range(64):
            u, v = x & 7, x >> 3
            phase = (u & v).bit_count() & 1
            phase ^= (z & 1) * ((u >> 0) & 1) * ((v >> 1) & 1)
            phase ^= ((z >> 1) & 1) * ((u >> 1) & 1) * ((v >> 2) & 1)
            phase ^= ((z >> 2) & 1) * ((u >> 2) & 1) * ((v >> 0) & 1)
            raw.append(1 - 2 * phase)
        image = walsh(raw)
        values = sorted(set(image))
        legal = all(abs(v) == 8 for v in image)
        assert legal == (z != 7)
        rows.append({'parameters': z, 'walsh_values': values, 'Boolean_slice': legal})
    return rows


parser = argparse.ArgumentParser()
parser.add_argument('input', type=Path)
parser.add_argument('--out', type=Path, required=True)
args = parser.parse_args()
raw = args.input.read_bytes()
data = json.loads(raw)
rows = [row for per_n in data['by_n'].values() for row in per_n['rows']]
nontrivial = sum(check_wrapper(row) for row in rows)
anchors = sum(check_existing_normal_forms(row) for row in rows)
result = {'status': 'passed', 'python': sys.version, 'platform': platform.platform(),
          'arithmetic': 'integer and F2', 'tolerance': 0, 'seed': None,
          'input_sha256': hashlib.sha256(raw).hexdigest(),
          'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
          'wrappers_checked': len(rows), 'one_heavy_wrappers': nontrivial,
          'existing_one_heavy_anchors_checked': anchors,
          'joint_counterexample': joint_counterexample()}
args.out.write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({k: v for k, v in result.items() if k not in ('joint_counterexample', 'python', 'platform')}))
