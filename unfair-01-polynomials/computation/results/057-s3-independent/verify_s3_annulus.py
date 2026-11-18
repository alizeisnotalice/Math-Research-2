#!/usr/bin/env python3
"""Verify S={3} surviving mothers with a 16th-degree integer norm."""

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import gzip
import hashlib
from itertools import groupby
import json
from math import gcd, isqrt
from pathlib import Path


def canonical(row):
    return json.dumps(row, separators=(",", ":")).encode()+b"\n"


def norm_polynomial(row):
    L, M, T, N, D, b, c, d, e = row
    doubled = (2, 2*L, 2*M, T, 2*b, 2*c, 2*d, 2*e, 2)
    p = [0]*17
    for i, value in enumerate(doubled):
        p[2*i] += value*value
        for j in range(i):
            p[i+j] += 2*value*doubled[j]
    p[6] -= D
    assert p[0] == p[-1] == 4 and all(x % 4 == 0 for x in p)
    return [x//4 for x in p]


def schur_descending(coefficients):
    """Return first failed degree, or zero for strict open-disk stability."""
    n = len(coefficients)-1
    a = coefficients
    while n:
        leading, constant = a[0], a[n]
        if abs(constant) >= abs(leading):
            return n
        reduction = [leading*a[i]-constant*a[n-i] for i in range(n)]
        content = gcd(*reduction)
        assert content > 0 and reduction[0] > 0
        a = [x//content for x in reduction]
        n -= 1
    return 0


def annulus(H):
    outer = schur_descending([H[i] << i for i in range(16, -1, -1)])
    if outer:
        return 16+outer
    inner = schur_descending([H[i] << (16-i) for i in range(17)])
    return 32+inner if inner else 0


def fixture_polynomial(factor):
    H = [0]*17
    for i, value in enumerate(factor):
        H[i] += value
        H[i+14] += value
    return H


def fixtures():
    assert annulus(fixture_polynomial((1, 0, 1))) == 0
    assert 17 <= annulus(fixture_polynomial((9, -15, 4))) <= 32  # roots 3, 3/4
    assert 17 <= annulus(fixture_polynomial((6, -11, 4))) <= 32  # roots 2, 3/4
    assert 33 <= annulus(fixture_polynomial((3, -11, 6))) <= 48  # roots 1/3, 3/2
    assert 33 <= annulus(fixture_polynomial((3, -8, 4))) <= 48   # roots 1/2, 3/2
    assert annulus(norm_polynomial((0,0,0,-6,24,0,0,0,0))) == 0
    assert annulus(norm_polynomial((0,0,0,-6,24,0,0,0,9))) != 0
    return 7


def verify_input_row(row):
    assert len(row) == 9 and all(type(x) is int for x in row), row
    L, M, T, N, D, b, c, d, e = row
    assert M == L*L and D == T*T-4*N and D > 0 and isqrt(D)**2 != D, row
    assert (-4 <= L <= 4 and -69 <= T <= 67 and -117 <= N <= 1189 and 5 <= D <= 468
            and -22 <= b <= 65 and -72 <= c <= 66 and -41 <= d <= 41 and -9 <= e <= 9), row


def read_tail_groups(manifest_root, entry, core):
    shard_path = manifest_root/entry["path"]
    payload = shard_path.read_bytes()
    assert (len(payload), hashlib.sha256(payload).hexdigest()) == (entry["bytes"],entry["sha256"])
    plain = gzip.decompress(payload)
    assert (len(plain), hashlib.sha256(plain).hexdigest()) == (entry["raw_bytes"],entry["raw_sha256"])
    groups = {}
    total = 0
    for line in plain.splitlines():
        values = json.loads(line)
        assert len(values) == 14 and tuple(values[:5]) == core
        key = tuple(values[5:7])
        assert key not in groups
        groups[key] = tuple(values[7:13])
        total += values[-1]
    assert len(groups) == entry["groups"] and total == entry["mothers"]
    return groups


def verify_core(job):
    rows, manifest_root, entry = job
    core = rows[0][:5]
    groups = read_tail_groups(manifest_root,entry,core) if entry is not None else None
    counts = Counter()
    failures = []
    trace = hashlib.sha256()
    for row in rows:
        if groups is not None:
            L = row[0]
            d, e = row[7:]
            dlo,dhi,elo,ehi,zlo,zhi = groups[row[5:7]]
            assert dlo <= d <= dhi and elo <= e <= ehi and zlo <= e-L*d <= zhi, row
        H = norm_polynomial(row)
        code = annulus(H)
        counts[code] += 1
        trace.update(canonical((row,H,code)))
        if code:
            failures.append((row,code))
    return core, counts, trace.hexdigest(), failures


def ingest(path):
    binary = path.read_bytes()
    plain = gzip.decompress(binary) if path.suffix == ".gz" else binary
    payload = hashlib.sha256()
    rows = []
    if plain:
        assert plain.endswith(b"\n")
    for line in plain.splitlines():
        row = tuple(json.loads(line))
        verify_input_row(row)
        payload.update(canonical(row))
        rows.append(row)
    return rows, dict(path=str(path), records=len(rows), bytes=len(binary),
                      sha256=hashlib.sha256(binary).hexdigest(), raw_bytes=len(plain),
                      raw_sha256=hashlib.sha256(plain).hexdigest(),
                      canonical_sha256=payload.hexdigest())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inputs", type=Path, nargs="+")
    ap.add_argument("--manifest", type=Path)
    ap.add_argument("--expected-count", type=int)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--allow-rejections", action="store_true")
    a = ap.parse_args()
    assert not a.out.exists()
    manifest = None
    entries = {}
    if a.manifest:
        raw = a.manifest.read_bytes()
        decoded = [json.loads(line) for line in raw.splitlines()]
        entries = {tuple(x["core"]):x for x in decoded}
        assert len(entries) == len(decoded)
        manifest = dict(path=str(a.manifest), bytes=len(raw), records=len(decoded),
                        sha256=hashlib.sha256(raw).hexdigest())
    n_fixtures = fixtures()
    rows, sources = [], []
    for path in a.inputs:
        partial, meta = ingest(path)
        rows.extend(partial)
        sources.append(meta)
    rows.sort()
    assert all(left < right for left,right in zip(rows,rows[1:])), "duplicate production result"
    if a.expected_count is not None:
        assert len(rows) == a.expected_count
    row_sha = hashlib.sha256()
    core_sha = hashlib.sha256()
    counts = Counter()
    failures = []
    checked_cores = 0
    root = a.manifest.parent if a.manifest else None
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        jobs = ((list(block), root, entries[core] if a.manifest else None)
                for core, block in groupby(rows, lambda row:row[:5]))
        for core, local, digest, rejected in pool.map(verify_core, jobs, chunksize=1):
            counts.update(local)
            core_sha.update(canonical((core,digest)))
            failures.extend(rejected)
            checked_cores += 1
    for row in rows:
        row_sha.update(canonical(row))
    assert sum(counts.values()) == len(rows)
    failure_blob = b"".join(map(canonical,failures))
    if failures:
        a.out.with_suffix(".rejections.jsonl").write_bytes(failure_blob)
    complete = not failures or a.allow_rejections
    summary = dict(status="verified" if complete else "mismatch", inputs=sources, rows=len(rows),
                   input_manifest=manifest, checked_cores=checked_cores, fixtures=n_fixtures,
                   failures=len(failures), failure_sha256=hashlib.sha256(failure_blob).hexdigest(),
                   annulus_codes=dict(sorted(counts.items())), canonical_rows_sha256=row_sha.hexdigest(),
                   canonical_norm_trace_sha256=core_sha.hexdigest(), workers=a.workers)
    a.out.write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    print(json.dumps(summary,sort_keys=True))
    if not complete:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
