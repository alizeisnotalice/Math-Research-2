#!/usr/bin/env python3
"""Construct explicit nonlinear switch partitions and small exact witnesses."""
import hashlib
import json
import platform
import time
from pathlib import Path


def basis(xs):
    pivots = {}
    answer = []
    for value in xs:
        x = value
        while x:
            p = x.bit_length() - 1
            if p in pivots:
                x ^= pivots[p]
            else:
                pivots[p] = x
                answer.append(value)
                break
    return answer


def fwht(values):
    a = list(values)
    step = 1
    while step < len(a):
        for i in range(0, len(a), 2 * step):
            for j in range(i, i + step):
                l, r = a[j], a[j + step]
                a[j], a[j + step] = l + r, l - r
        step *= 2
    return a


def construct(r):
    n = 6 + 2 * r
    size, q = 1 << n, 1 << (3 + r)
    blocks = [[] for _ in range(4 << r)]
    f0 = 0
    for x in range(size):
        z = [(x >> j) & 1 for j in range(6)]
        quad = (z[1] & z[2]) ^ (z[3] & z[4])
        a, b = z[0], z[5] ^ quad
        u, v = (x >> 6) & ((1 << r) - 1), x >> (6 + r)
        f = (z[1] & z[2]) ^ (a & (z[1] ^ z[3])) ^ (b & (z[2] ^ z[4])) ^ ((u & v).bit_count() & 1)
        f0 |= f << x
        blocks[a + 2 * b + 4 * v].append(x)
    assert all(len(a) == 2 * q for a in blocks)
    assert sorted(x for a in blocks for x in a) == list(range(size))
    generators = [(1 << size) - 1]
    generators += [sum(((x >> j) & 1) << x for x in range(size)) for j in range(n)]
    generators += [sum(1 << x for x in a) for a in blocks]
    bs = basis(generators)
    assert len(bs) == (4 << r) + r + 5
    # Each base partial Walsh support occupies a disjoint frequency set.
    union = set()
    partial = []
    for block in blocks:
        s = set(block)
        w = fwht([(1 - 2 * ((f0 >> x) & 1)) if x in s else 0 for x in range(size)])
        support = {i for i, c in enumerate(w) if c}
        assert len(support) == 2 * q and not support & union
        assert all(abs(w[i]) == q for i in support)
        union |= support
        partial.append(sorted(support))
    assert union == set(range(size))
    checked = 0
    if r == 0:
        # Complete first code; larger codes are independently checked remotely.
        for p in range(1 << len(bs)):
            phase = f0
            for j, b in enumerate(bs):
                if (p >> j) & 1:
                    phase ^= b
            assert all(abs(c) == q for c in fwht([1 - 2 * ((phase >> x) & 1) for x in range(size)]))
            checked += 1
    return {"r": r, "n": n, "q": q, "dimension": len(bs), "phase_count": 1 << len(bs),
            "f0_hex": format(f0, f'0{size // 4}x'),
            "basis_hex": [format(v, f'0{size // 4}x') for v in bs],
            "bits": "x0,x1,x2,x3,x4,x5,u[0:r],v[0:r]",
            "block_index": "a+2*b+4*v, a=x0, b=x5+x1*x2+x3*x4",
            "blocks": blocks, "base_partial_walsh_supports": partial,
            "local_full_phases_checked": checked}


def main():
    start = time.monotonic()
    root = Path(__file__).parent
    data = {"schema": "a28-r018-quadratic-graph-cases-v1", "cases": [construct(0), construct(1)],
            "python": platform.python_version(), "arithmetic": "exact integer/GF(2)",
            "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "elapsed_seconds": time.monotonic() - start}
    (root / "graph-cases.json").write_text(json.dumps(data, indent=2) + '\n')
    print(json.dumps({"dimensions": [c['dimension'] for c in data['cases']],
                      "phase_counts": [c['phase_count'] for c in data['cases']],
                      "local_full_phases": sum(c['local_full_phases_checked'] for c in data['cases']),
                      "elapsed_seconds": data['elapsed_seconds']}))


if __name__ == '__main__':
    main()
