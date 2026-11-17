#!/usr/bin/env python3
"""Independent portable exact certificate check, with all kernel Walsh slices."""
import argparse
import hashlib
import json
import platform
from pathlib import Path
import sys
import time


def parity(x):
    return bin(x).count('1') & 1


def rank(rows):
    b = {}
    for x in rows:
        while x:
            p = x.bit_length() - 1
            if p in b:
                x ^= b[p]
            else:
                b[p] = x
                break
    return len(b)


def fwht(a):
    a = a[:]
    step = 1
    while step < len(a):
        for start in range(0, len(a), 2 * step):
            for i in range(start, start + step):
                x, y = a[i], a[i + step]
                a[i], a[i + step] = x + y, x - y
        step *= 2
    return a


def verify(d):
    if d['seed_profile'] == 'base':
        q = [parity((x & 7) & (x >> 3)) for x in range(64)]
    elif d['seed_profile'] == 'mm-nonlinear':
        pi = (0, 1, 3, 2, 6, 7, 5, 4)
        q = [parity((x & 7) & pi[x >> 3]) ^ int((((x >> 3) & 3) == 3) ^ ((x >> 5) & 1)) for x in range(64)]
    else:
        raise AssertionError('unknown profile')
    assert all(abs(v) == 8 for v in fwht([1 - 2 * a for a in q]))
    all_rows = []
    assert len(d['matching_by_delta']) == 63
    for delta, matching in enumerate(d['matching_by_delta'], 1):
        assert len(matching) == 16
        positives = {(x, x ^ delta) for x in range(64) if x < (x ^ delta) and q[x] == q[x ^ delta]}
        negatives = {(x, x ^ delta) for x in range(64) if x < (x ^ delta) and q[x] != q[x ^ delta]}
        assert len(positives) == len(negatives) == 16
        seen_p, seen_n = set(), set()
        for a, b, c, e in matching:
            assert (a, b) in positives and (c, e) in negatives
            assert (a, b) not in seen_p and (c, e) not in seen_n
            seen_p.add((a, b)); seen_n.add((c, e))
            all_rows.append((1 << a) ^ (1 << b) ^ (1 << c) ^ (1 << e))
        assert seen_p == positives and seen_n == negatives
    assert all_rows == [int(x, 16) for x in d['incidence_rows_hex']]
    R = rank(all_rows)
    basis = [int(x, 16) for x in d['nullspace_basis_hex']]
    assert rank(basis) == len(basis) == 64 - R == d['best_dimV']
    assert R == d['best_rank_R']
    assert all(parity(b & row) == 0 for b in basis for row in all_rows)
    C = [(1 << 64) - 1] + [sum(((x >> j) & 1) << x for x in range(64)) for j in range(6)]
    assert rank(C) == 7 and rank(C + basis) == len(basis)
    assert rank([int(x, 16) for x in d['C_basis_hex']] + C) == 7
    columns = [sum(((b >> x) & 1) << j for j, b in enumerate(basis)) for x in range(64)]
    assert len(set(columns)) == 64
    assert rank(columns) == len(basis)
    signmask = sum(q[x] << x for x in range(64))
    assert all(parity(signmask & row) == 1 for row in all_rows)
    masks = [0]
    for b in basis:
        masks += [x ^ b for x in masks]
    for mask in masks:
        values = [1 - 2 * (q[x] ^ ((mask >> x) & 1)) for x in range(64)]
        assert all(abs(v) == 8 for v in fwht(values))
    return {'rows': len(all_rows), 'rank_R': R, 'dimV': len(basis),
            'affine_support_dimension': len(basis) - 1,
            'all_kernel_masks_checked': len(masks), 'C_checked_from_actual_coordinates': True,
            'frequency_columns_distinct': True, 'nonhomogeneous_syndrome_checked': True}


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('certificate', type=Path)
    p.add_argument('--out', type=Path, required=True); a = p.parse_args()
    begun = time.monotonic(); raw = a.certificate.read_bytes()
    result = verify(json.loads(raw))
    result.update(status='passed', python=sys.version, platform=platform.platform(),
                  elapsed_seconds=time.monotonic() - begun,
                  certificate_sha256=hashlib.sha256(raw).hexdigest(),
                  verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    a.out.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))
