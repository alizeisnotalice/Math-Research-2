#!/usr/bin/env python3
"""Exact CPU-only d=15 prefix plan; never imports a GPU worker."""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import os
import platform
import sys
import time
from collections import Counter
from pathlib import Path

PAIR = ((0, 0), (1, 0), (0, 1))


def digest(x: object) -> str:
    return hashlib.sha256(json.dumps(x, separators=(',', ':')).encode()).hexdigest()


def decode(code: str) -> tuple[list[int], list[int]]:
    a, b = [1], [1]
    for c in code:
        x, y = PAIR[ord(c) - ord('0')]
        a.append(x); b.append(y)
    return a, b


def convolution(a: list[int], b: list[int], end: int) -> list[int]:
    return [sum(a[i] * b[k-i] for i in range(max(0, k-len(b)+1), min(k, len(a)-1)+1)) for k in range(end+1)]


def decision(code: str, d: int) -> tuple[bool, tuple[str, ...]]:
    a, b = decode(code)
    if set(convolution(a, b, d-1)) - {0, 1} or sum(a[i] * b[d-i] for i in range(1, d)):
        return False, ()
    product = convolution(a, b, 2*d-2)
    reasons: list[str] = []
    if any(x > 1 for x in product): reasons.append('prefix_product_coefficient_gt_1')
    # At i=d-1 the formal baseline has degree 2d-1, beyond the product of
    # two (d-1)-prefixes; its coefficient is exactly zero.
    if any((product[d+i] if d+i < len(product) else 0) >= 1 and a[i] + b[i] > 0 for i in range(1, d)):
        reasons.append('forced_positive_fractional_cross_at_d_plus_i')
    return True, tuple(reasons)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument('--out', type=Path, required=True); p.add_argument('--d', type=int, default=15)
    p.add_argument('--m', type=int, default=64); p.add_argument('--n', type=int, default=104)
    p.add_argument('--minimum-starts', type=int, default=8192); p.add_argument('--seed', type=int, default=2026111701)
    p.add_argument('--margin', type=float, default=.001); a = p.parse_args()
    assert Path('/Volumes/data').is_mount() and a.d == 15 and a.m > a.d and a.n > a.d
    assert a.minimum_starts >= 8192 and a.minimum_starts % 8192 == 0 and 0 < a.margin < .5
    began = time.monotonic(); valid: list[str] = []; survivors: list[str] = []; rejected = Counter()
    for chars in itertools.product('012', repeat=a.d-1):
        code = ''.join(chars); legal, reasons = decision(code, a.d)
        if not legal: continue
        valid.append(code)
        if reasons: rejected.update(reasons)
        else: survivors.append(code)
    assert survivors
    starts = max(a.minimum_starts, ((len(survivors)+8191)//8192)*8192)
    assigned = [survivors[i % len(survivors)] for i in range(starts)]
    counts = Counter(assigned)
    assert set(counts) == set(survivors) and min(counts.values()) >= 1
    out = {'schema':'d15-necessary-prefix-plan-v1', 'status':'cpu_generated_gpu_not_started',
           'parameters':{'d':a.d,'m':a.m,'n':a.n,'minimum_starts':a.minimum_starts,'coverage_plan_starts':starts,
                         'margin':a.margin,'seed':a.seed,'suggested_steps':12000,'suggested_restart_every':2000},
           'raw_ternary_count':3**(a.d-1), 'valid_low_prefix_count':len(valid), 'survivor_count':len(survivors),
           'valid_low_prefix_codes_sha256':digest(valid),'survivor_codes_sha256':digest(survivors),
           'strict_filter_rejected_by_reason_counts':dict(rejected), 'survivor_codes':survivors,
           'coverage_assignment_codes':assigned,
           'per_prefix_assignment_counts':[{'code':x,'count':counts[x]} for x in survivors],
           'elapsed_seconds':time.monotonic()-began,'interpreter':sys.version,'platform':platform.platform(),
           'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           'claim_boundary':'Exact finite low-prefix/necessary-filter plan only. It neither starts GPU nor proves nonexistence or factorization.'}
    a.out.parent.mkdir(parents=True, exist_ok=True); a.out.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:out[k] for k in ('raw_ternary_count','valid_low_prefix_count','survivor_count','elapsed_seconds')},sort_keys=True))


if __name__ == '__main__':
    main()
