#!/usr/bin/env python3
"""Independent semantic audit of boundary-bitset rank traces.

This reader does not call Metal or reuse the Swift CPU reference.  It checks
every emitted GPU state, allowed-count vector, CDF, keyed parent choice, and
mutation against Python integers and a direct grid-neighbour computation.
"""
import json
import sys
from pathlib import Path

MASK = (1 << 64) - 1
STEP = 0x9E3779B97F4A7C15
C1 = 0xBF58476D1CE4E5B9
C2 = 0x94D049BB133111EB
K1 = 0x9E3779B97F4A7C15
K2 = 0xD1B54A32D192ED03


def mix(x):
    x = ((x ^ (x >> 30)) * C1) & MASK
    x = ((x ^ (x >> 27)) * C2) & MASK
    return x ^ (x >> 31)


def below(state, bound):
    threshold = ((-bound) & MASK) % bound
    while True:
        state = (state + STEP) & MASK
        value = mix(state)
        if value >= threshold:
            return value % bound


def key(seed, rank, child, domain):
    return seed ^ (rank << 32) ^ ((child * K1) & MASK) ^ ((domain * K2) & MASK)


def split(words, particle):
    base = particle * 66
    return words[base:base + 32], words[base + 32:base + 64], words[base + 64:base + 66]


def occupied(words):
    return {32 * w + j for w, word in enumerate(words) for j in range(32) if word >> j & 1}


def boundary(n, occ):
    ans = set()
    for v in occ:
        x, y = v % n, v // n
        for q in (v - 1 if x else -1, v + 1 if x + 1 < n else -1,
                  v - n if y else -1, v + n if y + 1 < n else -1):
            if q >= 0 and q not in occ:
                ans.add(q)
    return ans


def mask_from_set(items):
    out = [0] * 32
    for v in items:
        out[v >> 5] |= 1 << (v & 31)
    return out


def candidates(n, occ, bd, tau, rank):
    available = set(range(n * n)) - occ
    if rank == 1:
        return sorted(available)
    return sorted(available - bd if rank == tau else bd)


def verify(path):
    rows = [json.loads(line) for line in path.read_text().splitlines()]
    meta = rows[0]
    assert meta['schema'] == 'bitset-rank-trace-v1'
    n, tau, p = meta['n'], meta['tau'], meta['particles']
    seed = int(meta['seed_hex'], 16)
    ranks = rows[1:-1]
    assert rows[-1] == {'kind': 'audit_end', 'success': True}
    assert len(ranks) == n * n
    previous = [[0] * 66 for _ in range(p)]
    checks = 0
    dead = False
    for row in ranks:
        rank = row['rank']
        assert row['matched'] and row['cpu_boundary_recomputed']
        pre = previous
        if dead:
            assert row['allowed'] == [0] * p
            assert row['groups'] == [0] * ((p + 255) // 256)
            assert row['cdf'] == [0] * ((p + 255) // 256)
            assert row['total'] == 0
            if rank == n * n:
                assert row['parents'] is None
            else:
                assert row['parents'] == [0xffffffff] * p
            assert row['state_words'] == [word for state in pre for word in state]
            checks += p
            continue
        expected_allowed = []
        allowed_lists = []
        for state in pre:
            occ_words, bd_words, _ = split(state, 0)
            occ = occupied(occ_words)
            bd = boundary(n, occ)
            assert mask_from_set(bd) == bd_words
            assert not (set(range(n * n, 1024)) & occ)
            assert not (set(range(n * n, 1024)) & occupied(bd_words))
            choices = candidates(n, occ, bd, tau, rank)
            allowed_lists.append(choices)
            expected_allowed.append(len(choices))
        assert row['allowed'] == expected_allowed
        groups = [sum(expected_allowed[i:i + 256]) for i in range(0, p, 256)]
        assert row['groups'] == groups
        cdf, total = [], 0
        for g in groups:
            total += g
            cdf.append(total)
        assert row['cdf'] == cdf and row['total'] == total
        if total == 0:
            if rank == n * n:
                assert row['parents'] is None
            else:
                assert row['parents'] == [0xffffffff] * p
            assert row['state_words'] == [word for state in pre for word in state]
            dead = True
            checks += p
            continue
        if rank == n * n:
            assert row['parents'] is None
            for state in (row['state_words'][66 * i:66 * (i + 1)] for i in range(p)):
                occ_words, bd_words, _ = split(state, 0)
                assert len(occupied(occ_words)) == n * n - 1
                assert mask_from_set(boundary(n, occupied(occ_words))) == bd_words
            checks += p
            continue
        parents = row['parents']
        words = row['state_words']
        assert len(parents) == p and len(words) == 66 * p
        for child, parent in enumerate(parents):
            first = below(key(seed, rank, child, 1), total)
            u = (child * total + first) // p if meta['method'] == 'strat' else first
            cumulative = 0
            expected_parent = None
            for i, weight in enumerate(expected_allowed):
                if u < cumulative + weight:
                    expected_parent = i
                    break
                cumulative += weight
            assert parent == expected_parent
            old_o, old_b, old_roots = split(pre[parent], 0)
            pick = below(key(seed, rank, child, 2), expected_allowed[parent])
            new_vertex = allowed_lists[parent][pick]
            new = words[66 * child:66 * (child + 1)]
            new_o, new_b, new_roots = split(new, 0)
            old_occ, new_occ = occupied(old_o), occupied(new_o)
            assert new_occ == old_occ | {new_vertex}
            assert mask_from_set(boundary(n, new_occ)) == new_b
            roots = list(old_roots)
            if rank == 1:
                roots[0] = new_vertex
            if rank == tau:
                roots[1] = new_vertex
            assert new_roots == roots
            checks += 1
        previous = [words[66 * i:66 * (i + 1)] for i in range(p)]
    return {'file': path.name, 'n': n, 'tau': tau, 'particles': p,
            'ranks': len(ranks), 'state_parent_rng_checks': checks}


root = Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).parent
targets = (
    [root / f'trace-n{n}-tau2-p17-fixed.jsonl' for n in (12, 16, 24, 32)]
    + sorted(root.glob('trace-final-n*-tau*-p17.jsonl'))
)
result = {'status': 'passed', 'method': 'independent Python state, boundary, candidate, keyed-RNG and parent audit',
          'rows': [verify(path) for path in targets]}
(root / 'independent-trace-audit.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
