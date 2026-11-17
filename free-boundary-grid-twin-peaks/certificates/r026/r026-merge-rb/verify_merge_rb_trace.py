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
from fractions import Fraction

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


def groups(values: list[int], width: int = 256) -> list[int]:
    return [sum(values[i:i + width]) for i in range(0, len(values), width)]


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit('usage: verify_rankbiased_trace.py raw.jsonl trace.json')
    raw = [json.loads(x) for x in pathlib.Path(sys.argv[1]).read_text().splitlines() if x.strip()]
    metadata = raw[0]
    assert metadata['implementation'] == 'r026-tau2-firstmerge-rb-terminal'
    assert metadata['base_rankbiased_source_sha256'] == '60746225d734d03a5567c0b2e45581e89cd8a39599605f6e1a4397b612e23334'
    assert metadata['tau_lo'] == metadata['tau_hi'] == 2
    trace = json.loads(pathlib.Path(sys.argv[2]).read_text())
    replicate = next(x for x in raw if x['kind'] == 'replicate')
    n, P, tau, L = trace['n'], trace['P'], trace['tau'], trace['L']
    # The trace's outer header remains copied r025 provenance; terminal r026
    # payload and raw r026 metadata are checked below and are authoritative.
    assert P == metadata['particles'] and tau == 2 and len(trace['rows']) == n * n
    seed = int(trace['seed_hex'], 16)
    words = [tuple() for _ in range(P)]
    moves: dict[str, int] = {}
    audited = []
    for rank, row in enumerate(trace['rows'], 1):
        assert row['rank'] == rank
        weights = [len(allowed(n, tau, w)) for w in words]
        assert row['allowed'] == weights
        running = sum(weights)
        cdf = []
        subtotal = 0
        for i in range(0, P, 256):
            subtotal += sum(weights[i:i + 256])
            cdf.append(subtotal)
        assert row['cdf'] == [str(x) for x in cdf]
        if row.get('terminal'):
            assert rank == n * n
            assert row['prefix'] == [list(w) for w in words]
            assert row['state'] == [state(n, w) for w in words]
            # The terminal trace intentionally stores N-1 entries.  Append
            # the unique omitted point only for this independent audit.
            Lsum = LDsum = LFsum = LGsum = L2sum = 0
            msum = mdsum = mfsum = mgsum = Fraction(0)
            mhist = [Fraction(0) for _ in range(2 * n - 1)]
            hist = [0] * (2 * n - 1)
            Lvals: list[int] = []; LDvals: list[int] = []; LFvals: list[int] = []; LGvals: list[int] = []; L2vals: list[int] = []
            d2sum = fsum = gsum = 0
            Dvals: list[int] = []; Fvals: list[int] = []; Gvals: list[int] = []
            payload = row['first_merge_rb_payload_ab_r1_r2']
            assert row['first_merge_rb_bad'] == 0 and len(payload) == 4 * P
            for particle, word in enumerate(words):
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
                Lvals.append(ell); LDvals.append(ell * d2); LFvals.append(ell * f0); LGvals.append(ell * g0); L2vals.append(ell * ell)
                d2sum += d2; fsum += f0; gsum += g0
                Dvals.append(d2); Fvals.append(f0); Gvals.append(g0)
                # Independent complete-word first-merge scan: components are
                # represented simply by the two growing sets up to first join.
                A, B = {vx}, {vy}
                for v in full[2:]:
                    ta = any(v in neighbours(n, u) for u in A)
                    tb = any(v in neighbours(n, u) for u in B)
                    if ta and tb:
                        break
                    assert ta != tb
                    (A if ta else B).add(v)
                else:
                    raise AssertionError('first merge missing')
                a, b = len(A), len(B)
                assert payload[4 * particle:4 * particle + 4] == [a, b, vx, vy]
                # Exact (not floating) per-payload cone check.  Its proof is
                # q>=d and 4ab<=q^2, so it is structural rather than an
                # accuracy diagnostic.
                q = a + b
                assert q * (q - 1) * md >= 4 * a * b * (md - 1)
                m = Fraction((a + b) * (a + b - 1), 2 * a * b)
                msum += m; mdsum += m * d2; mfsum += m * f0; mgsum += m * g0; mhist[md] += m
            assert row['tau2_lift_bad'] == row['tau2_distance_bad'] == 0
            assert groups(Lvals) == [int(x) for x in row['tau2_lift_L_groups']]
            assert groups(LDvals) == [int(x) for x in row['tau2_lift_LD_groups']]
            assert groups(LFvals) == [int(x) for x in row['tau2_lift_LF_groups']]
            assert groups(LGvals) == [int(x) for x in row['tau2_lift_LG_groups']]
            assert [groups([ell if abs(word[0] % n - word[1] % n) + abs(word[0] // n - word[1] // n) == d else 0 for word, ell in zip(words, Lvals)]) for d in range(2 * n - 1)] == list(map(list, zip(*[[int(x) for x in group] for group in row['tau2_distance_L_hist_groups']])))
            assert groups(L2vals) == [int(x) for x in row['tau2_distance_L2_groups']]
            assert LFsum - 2 * LGsum == 4 * LDsum
            assert groups(Dvals) == [int(x) for x in row['moment_groups']]
            assert groups(Fvals) == [int(x) for x in row['root_F_groups']]
            assert groups(Gvals) == [int(x) for x in row['root_G_groups']]
            for field, value in (
                ('tau2_lift_L_groups', groups(Lvals)), ('tau2_lift_L_distance2_groups', groups(LDvals)),
                ('tau2_lift_L_F0_groups', groups(LFvals)), ('tau2_lift_L_G0_groups', groups(LGvals)),
                ('tau2_distance_L2_groups', groups(L2vals)),
            ):
                assert replicate[field] == value
            assert abs(float(replicate['U_over_Z']) - d2sum / P / (n * n)) < 1e-14
            assert abs(float(replicate['F_over_Z']) - fsum / P) < 1e-14
            assert abs(float(replicate['G_over_Z']) - gsum / P) < 1e-14
            assert abs(float(replicate['tau2_lift_scaled_mean_L']) - Lsum / P) < 1e-14
            assert abs(float(replicate['tau2_lift_scaled_mean_L_distance2']) - LDsum / P) < 1e-14
            assert abs(float(replicate['tau2_lift_scaled_mean_L_F0']) - LFsum / P) < 1e-14
            assert abs(float(replicate['tau2_lift_scaled_mean_L_G0']) - LGsum / P) < 1e-14
            assert replicate['tau2_lift_group_LF_minus_2LG_minus_4LD'] == [0] * len(groups(Lvals))
            assert replicate['tau2_distance_scaled_mean_L_hist'] == [h / P for h in hist]
            assert replicate['tau2_distance_scaled_mean_L2'] == L2sum / P
            assert mfsum - 2 * mgsum == 4 * mdsum
            left, right = Fraction(P, 1) / msum, sum(mhist[d] * Fraction(d, 2 * (d - 1)) for d in range(2, 2 * n - 1)) / msum
            assert left <= right
            audited.append({'rank': rank, 'terminal': True, 'parents': None})
            continue
        expected_parents, after = [], []
        for particle in range(P):
            z = key(seed, rank, particle, 1)
            z, u = below(z, running)
            u = (particle * running + u) // P if trace['method'] == 'strat' else u
            group = next(i for i, cut in enumerate(cdf) if u < cut)
            low, prior = 256 * group, (cdf[group - 1] if group else 0)
            subtotal = prior
            parent = low
            for candidate in range(low, min(low + 256, P)):
                if u < subtotal + weights[candidate]:
                    parent = candidate
                    break
                subtotal += weights[candidate]
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
        'first_merge_rb': 'per-particle GPU payload matches independent full-word A/B scan including terminal omitted vertex; Fraction mH identity and finite moment inequality checked',
        'proposal_categories': moves,
        'terminal_full_word_reconstruction': 'passed: appended unique omitted vertex and verified K=2, L, d/L fields, root fields, and all group/raw identities',
        'caveat': 'Trace equality validates this finite implementation path only; it is not a mixing or SMC-accuracy certificate.'
    }
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
