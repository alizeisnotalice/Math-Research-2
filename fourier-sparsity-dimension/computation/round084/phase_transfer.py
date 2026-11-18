"""Exact support and affine-rank audit for unconstrained phase functions."""
import json
import multiprocessing
import platform
import random
import sys


def walsh(values):
    w = values.copy()
    step = 1
    while step < len(w):
        for j in range(0, len(w), 2 * step):
            for i in range(j, j + step):
                a, b = w[i], w[i + step]
                w[i], w[i + step] = a + b, a - b
        step *= 2
    return w


def basis(values):
    pivots = {}
    result = []
    for value in values:
        v = value
        while v and v.bit_length() in pivots:
            v ^= pivots[v.bit_length()]
        if v:
            pivots[v.bit_length()] = v
            result.append(value)
    return result


def batch(seed):
    rng = random.Random(seed)
    sparse = 0
    cancellations = 0
    for case in range(96):
        g = [rng.choice([-1, 1]) for _ in range(16)]
        if case % 3 == 0:
            g = [(-1 if (x >> (1 + (x & 1))) & 1 else 1) for x in range(16)]
        shifts = [rng.randrange(16) for _ in range(8)]
        signs = [rng.choice([-1, 1]) for _ in shifts]
        parent = walsh([signs[u] * g[x ^ shifts[u]] for u in range(8) for x in range(16)])
        gw = walsh(g)
        support = [e for e in range(16) if gw[e]]
        sparse += len(support) < 16
        rng.shuffle(support)
        chosen = {x & 15 for x in basis([16 + e for e in support])}
        candidate = []
        for e in support:
            p = walsh([s * (-1)**((t & e).bit_count() % 2) for t, s in zip(shifts, signs)])
            active = [128 + (u << 4) + e for u in range(8) if p[u]]
            assert all(parent[(u << 4) + e] == p[u] * gw[e] for u in range(8))
            if e in chosen:
                candidate.extend(basis(active))
            else:
                cancellations += 8 - len(active)
        whole = [128 + j for j, a in enumerate(parent) if a]
        assert len(basis(candidate)) == len(basis(whole))
    return {'seed': seed, 'cases': 96, 'sparse': sparse, 'absent_in_dependent_bands': cancellations}


if __name__ == '__main__':
    jobs, seed = map(int, sys.argv[1:])
    with multiprocessing.get_context('fork').Pool(jobs) as pool:
        result = pool.map(batch, range(seed, seed + jobs))
    print(json.dumps(dict(python=platform.python_version(), jobs=jobs, results=result), indent=2))
