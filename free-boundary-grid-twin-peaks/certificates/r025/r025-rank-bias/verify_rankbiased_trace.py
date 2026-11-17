#!/usr/bin/env python3
"""Replay the small r025 GPU trace with integer SplitMix64 arithmetic.

This verifier is intentionally independent of Metal state storage.  It
reconstructs the fully adapted resampling, rank-biased fixed-prefix move, and
the occupancy/boundary state from the trace seed.  It is only sized for the
small audit trace; it is not an SMC production implementation.
"""
from __future__ import annotations

import json
import pathlib
import sys

MASK = (1 << 64) - 1
GOLD = 0x9E3779B97F4A7C15
M1 = 0xBF58476D1CE4E5B9
M2 = 0x94D049BB133111EB
C1 = 0x9E3779B97F4A7C15
C2 = 0xD1B54A32D192ED03


def nxt(x: int) -> tuple[int, int]:
    x = (x + GOLD) & MASK
    z = x
    z = ((z ^ (z >> 30)) * M1) & MASK
    z = ((z ^ (z >> 27)) * M2) & MASK
    return x, (z ^ (z >> 31)) & MASK


def below(x: int, n: int) -> tuple[int, int]:
    threshold = ((-n) & MASK) % n
    while True:
        x, z = nxt(x)
        if z >= threshold:
            return x, z % n


def key(seed: int, rank: int, particle: int, draw: int) -> int:
    return (seed ^ (rank << 32) ^ ((particle * C1) & MASK) ^ ((draw * C2) & MASK)) & MASK


def neighbours(n: int, v: int):
    x, y = v % n, v // n
    if x:
        yield v - 1
    if x + 1 < n:
        yield v + 1
    if y:
        yield v - n
    if y + 1 < n:
        yield v + n


def boundary(n: int, word: tuple[int, ...]) -> set[int]:
    seen = set(word)
    return {w for v in seen for w in neighbours(n, v) if w not in seen}


def allowed(n: int, tau: int, word: tuple[int, ...]) -> list[int]:
    seen = set(word)
    rank = len(word) + 1
    outside = [v for v in range(n * n) if v not in seen]
    if rank == 1:
        return outside
    b = boundary(n, word)
    if rank == tau:
        return [v for v in outside if v not in b]
    return [v for v in outside if v in b]


def valid(n: int, tau: int, word: tuple[int, ...]) -> bool:
    if len(word) != len(set(word)):
        return False
    seen: set[int] = set()
    births = []
    for rank, v in enumerate(word, 1):
        if not any(w in seen for w in neighbours(n, v)):
            births.append(rank)
        seen.add(v)
    return tuple(births) == ((1,) if len(word) < tau else (1, tau))


def birth_ranks(n: int, word: tuple[int, ...]) -> tuple[int, ...]:
    seen: set[int] = set()
    out = []
    for rank, v in enumerate(word, 1):
        if not any(w in seen for w in neighbours(n, v)):
            out.append(rank)
        seen.add(v)
    return tuple(out)


def move(n: int, tau: int, word: tuple[int, ...], rank: int, particle: int, seed: int, steps: int, moves: dict[str, int]) -> tuple[int, ...]:
    out = word
    for step in range(steps):
        z = key(seed, rank, particle, step + 3)
        i = 0
        if rank >= 2:
            z, branch = below(z, 2)
            if branch:
                z, i = below(z, rank)
            else:
                z, i = below(z, 2)
        z, v = below(z, n * n)
        y = list(out)
        try:
            j = y.index(v)
            y[i], y[j] = y[j], y[i]
            if i == j:
                category = 'self_active'
            elif (i < min(2, rank)) != (j < min(2, rank)):
                category = 'root_nonroot_active'
            elif i < min(2, rank):
                category = 'root_root_active'
            else:
                category = 'nonroot_active'
        except ValueError:
            y[i] = v
            category = 'root_outside' if i < min(2, rank) else 'nonroot_outside'
        candidate = tuple(y)
        moves[category + '_proposed'] = moves.get(category + '_proposed', 0) + 1
        if valid(n, tau, candidate):
            out = candidate
            moves[category + '_accepted'] = moves.get(category + '_accepted', 0) + 1
        else:
            moves[category + '_rejected'] = moves.get(category + '_rejected', 0) + 1
    return out


def state(n: int, word: tuple[int, ...]) -> list[object]:
    occ = [0] * 32
    ext = [0] * 32
    seen = set(word)
    for v in seen:
        occ[v >> 5] |= 1 << (v & 31)
    for v in boundary(n, word):
        ext[v >> 5] |= 1 << (v & 31)
    return [occ, ext, word[0] if word else 0, word[1] if len(word) >= 2 else 0]


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit('usage: verify_rankbiased_trace.py raw.jsonl trace.json')
    raw = [json.loads(x) for x in pathlib.Path(sys.argv[1]).read_text().splitlines() if x.strip()]
    metadata = raw[0]
    assert metadata['implementation'] == 'r025-tau2-rankbiased-prefix-mh'
    assert metadata['mh_move_kernel'] == 'r>=2: half uniform {0,1}, half uniform [0,r); r=1: rank 0'
    assert metadata['tau_lo'] == metadata['tau_hi'] == 2
    trace = json.loads(pathlib.Path(sys.argv[2]).read_text())
    replicate = next(x for x in raw if x['kind'] == 'replicate')
    n, P, tau, L = trace['n'], trace['P'], trace['tau'], trace['L']
    assert trace['implementation'] == metadata['implementation']
    assert P == metadata['particles'] and tau == 2 and len(trace['rows']) == n * n
    seed = int(trace['seed_hex'], 16)
    words = [tuple() for _ in range(P)]
    moves: dict[str, int] = {}
    audited = []
    for rank, row in enumerate(trace['rows'], 1):
        assert row['rank'] == rank
        weights = [len(allowed(n, tau, w)) for w in words]
        assert row['allowed'] == weights
        cdf = []
        running = 0
        for weight in weights:
            running += weight
            cdf.append(running)
        assert row['cdf'] == [str(running)]
        if row.get('terminal'):
            assert rank == n * n
            assert row['prefix'] == [list(w) for w in words]
            assert row['state'] == [state(n, w) for w in words]
            # The terminal trace intentionally stores N-1 entries.  Append
            # the unique omitted point only for this independent audit.
            Lsum = LDsum = LFsum = LGsum = L2sum = 0
            hist = [0] * (2 * n - 1)
            d2sum = fsum = gsum = 0
            for word in words:
                missing = next(v for v in range(n * n) if v not in word)
                full = word + (missing,)
                assert birth_ranks(n, full) == (1, 2)
                vx, vy = full[:2]
                assert abs(vx % n - vy % n) + abs(vx // n - vy // n) >= 2
                best = min(k + 1 for k, w in enumerate(word)
                           if abs(vy % n - w % n) + abs(vy // n - w // n) == 1)
                assert 2 < best < n * n
                ell = best - 2
                dx, dy = vx % n - vy % n, vx // n - vy // n
                d2 = dx * dx + dy * dy
                md = abs(dx) + abs(dy)
                ax, ay = 2 * (vx % n) - (n - 1), 2 * (vx // n) - (n - 1)
                bx, by = 2 * (vy % n) - (n - 1), 2 * (vy // n) - (n - 1)
                f0, g0 = ax * ax + ay * ay + bx * bx + by * by, ax * bx + ay * by
                Lsum += ell; LDsum += ell * d2; LFsum += ell * f0; LGsum += ell * g0
                L2sum += ell * ell; hist[md] += ell
                d2sum += d2; fsum += f0; gsum += g0
            assert row['tau2_lift_bad'] == row['tau2_distance_bad'] == 0
            assert [Lsum] == [int(x) for x in row['tau2_lift_L_groups']]
            assert [LDsum] == [int(x) for x in row['tau2_lift_LD_groups']]
            assert [LFsum] == [int(x) for x in row['tau2_lift_LF_groups']]
            assert [LGsum] == [int(x) for x in row['tau2_lift_LG_groups']]
            assert [hist] == [[int(x) for x in group] for group in row['tau2_distance_L_hist_groups']]
            assert [L2sum] == [int(x) for x in row['tau2_distance_L2_groups']]
            assert LFsum - 2 * LGsum == 4 * LDsum
            assert [d2sum] == [int(x) for x in row['moment_groups']]
            assert [fsum] == [int(x) for x in row['root_F_groups']]
            assert [gsum] == [int(x) for x in row['root_G_groups']]
            for field, value in (
                ('tau2_lift_L_groups', [Lsum]), ('tau2_lift_L_distance2_groups', [LDsum]),
                ('tau2_lift_L_F0_groups', [LFsum]), ('tau2_lift_L_G0_groups', [LGsum]),
                ('tau2_distance_L_hist_groups', [hist]), ('tau2_distance_L2_groups', [L2sum]),
            ):
                assert replicate[field] == value
            assert abs(float(replicate['U_over_Z']) - d2sum / P / (n * n)) < 1e-14
            assert abs(float(replicate['F_over_Z']) - fsum / P) < 1e-14
            assert abs(float(replicate['G_over_Z']) - gsum / P) < 1e-14
            assert abs(float(replicate['tau2_lift_scaled_mean_L']) - Lsum / P) < 1e-14
            assert abs(float(replicate['tau2_lift_scaled_mean_L_distance2']) - LDsum / P) < 1e-14
            assert abs(float(replicate['tau2_lift_scaled_mean_L_F0']) - LFsum / P) < 1e-14
            assert abs(float(replicate['tau2_lift_scaled_mean_L_G0']) - LGsum / P) < 1e-14
            assert replicate['tau2_lift_group_LF_minus_2LG_minus_4LD'] == [0]
            assert replicate['tau2_distance_scaled_mean_L_hist'] == [h / P for h in hist]
            assert replicate['tau2_distance_scaled_mean_L2'] == L2sum / P
            audited.append({'rank': rank, 'terminal': True, 'parents': None})
            continue
        expected_parents, after = [], []
        for particle in range(P):
            z = key(seed, rank, particle, 1)
            z, u = below(z, running)
            u = (particle * running + u) // P if trace['method'] == 'strat' else u
            parent = next(i for i, cut in enumerate(cdf) if u < cut)
            expected_parents.append(parent)
            z = key(seed, rank, particle, 2)
            z, k = below(z, weights[parent])
            proposal = words[parent] + (allowed(n, tau, words[parent])[k],)
            after.append(move(n, tau, proposal, rank, particle, seed, L, moves))
        assert row['parents'] == expected_parents
        assert row['prefix'] == [list(w) for w in after]
        assert row['state'] == [state(n, w) for w in after]
        words = after
        audited.append({'rank': rank, 'terminal': False, 'parents': expected_parents})
    out = {
        'passed': True,
        'scope': 'Exact Python replay of resample, append, rank-biased move, and state fields for one GPU trace',
        'n': n, 'P': P, 'tau': tau, 'move_steps': L,
        'ranks_checked': len(audited),
        'rank_selection': metadata['mh_rank_weights'],
        'proposal_categories': moves,
        'terminal_full_word_reconstruction': 'passed: appended unique omitted vertex and verified K=2, L, d/L fields, root fields, and all group/raw identities',
        'caveat': 'Trace equality validates this finite implementation path only; it is not a mixing or SMC-accuracy certificate.'
    }
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
