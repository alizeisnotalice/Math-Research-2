#!/usr/bin/env python3
"""MPS-assisted exact search for common sign codes of 6-variable MM bent spectra.

The only floating-point operation is a 64 by 64 Walsh matrix product on MPS.
Every reported candidate is rechecked using an integer Walsh transform with no
tolerance.  The search family consists of signs obtained by flipping an
indicator of a three-dimensional affine flat in F_2^6.
"""
import argparse
import hashlib
import itertools
import json
import os
from pathlib import Path
import platform
import random
import sys
import time

os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "0")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")

import numpy as np
import torch

N = 6
POINTS = 1 << N
TARGET = 1 << (N // 2)


def parity(x):
    return x.bit_count() & 1


def rank(vectors):
    pivots = {}
    for vector in vectors:
        x = vector
        while x:
            p = x.bit_length() - 1
            if p in pivots:
                x ^= pivots[p]
            else:
                pivots[p] = x
                break
    return len(pivots)


def rref(vectors):
    pivots = {}
    for vector in vectors:
        x = vector
        while x:
            p = x.bit_length() - 1
            if p in pivots:
                x ^= pivots[p]
            else:
                pivots[p] = x
                break
    return pivots


def reduce_by(x, pivots):
    for p in sorted(pivots, reverse=True):
        if (x >> p) & 1:
            x ^= pivots[p]
    return x


def fwht(values):
    values = list(values)
    step = 1
    while step < len(values):
        for start in range(0, len(values), 2 * step):
            for j in range(start, start + step):
                a, b = values[j], values[j + step]
                values[j], values[j + step] = a + b, a - b
        step *= 2
    return values


def walsh_matrix():
    return np.array([[1 if not parity(x & y) else -1 for y in range(POINTS)]
                     for x in range(POINTS)], dtype=np.int16)


def enumerate_affine_3flats():
    """Return one 64-bit indicator for every 3-flat, exactly 11160 masks."""
    subspaces = set()
    nonzero = range(1, POINTS)
    for a in nonzero:
        for b in range(a + 1, POINTS):
            if b == a:
                continue
            for c in range(b + 1, POINTS):
                if rank([a, b, c]) != 3:
                    continue
                basis = rref([a, b, c])
                elems = {0}
                for v in basis.values():
                    elems |= {z ^ v for z in list(elems)}
                subspaces.add(tuple(sorted(elems)))
    masks = set()
    for subspace in subspaces:
        rep = min(range(POINTS), key=lambda x: tuple(sorted(x ^ y for y in subspace)))
        # The preceding rep is a canonical coset only when it is minimum in it.
        # Enumerate all cosets, then set-deduplicate across equivalent reps.
        for translate in range(POINTS):
            points = tuple(sorted(translate ^ y for y in subspace))
            if points[0] != translate and translate != min(points):
                continue
            masks.add(sum(1 << x for x in points))
    masks = sorted(masks)
    if len(subspaces) != 1395 or len(masks) != 11160:
        raise AssertionError((len(subspaces), len(masks)))
    bits = np.array([[(m >> x) & 1 for x in range(POINTS)] for m in masks], dtype=np.int16)
    return masks, bits


def mm_signs(permutation, h):
    """Raw (= eight times normalized) 6-variable Maiorana--McFarland signs."""
    return np.array([1 if not parity((x & 7) & permutation[x >> 3]) ^ ((h >> (x >> 3)) & 1)
                     else -1 for x in range(POINTS)], dtype=np.int16)


def c_space():
    generators = [(1 << POINTS) - 1]
    generators += [sum(((x >> j) & 1) << x for x in range(POINTS)) for j in range(N)]
    pivots = rref(generators)
    if len(pivots) != N + 1:
        raise AssertionError("C has unexpected rank")
    return pivots


def candidate_masks(base, flat_bits, hmat):
    signs = base[None, :] * (1 - 2 * flat_bits)
    raw = signs @ hmat
    legal = np.all(np.abs(raw) == TARGET, axis=1)
    return signs, raw, legal


def mps_candidates(signs, h_mps):
    # Includes upload and download in the measured interval.  h_mps is resident.
    began = time.perf_counter()
    got = torch.from_numpy(signs.astype(np.float32, copy=False)).to("mps") @ h_mps
    torch.mps.synchronize()
    raw = got.cpu().numpy()
    elapsed = time.perf_counter() - began
    return raw, elapsed


def pairwise_trap(hmat, h_mps):
    """A three-mask test where every pair but not the full span is legal.

    q_z(u,v)=u.(I+z1E01+z2E12+z3E20)v has determinant
    1+z1*z2*z3 over F_2.  Thus 111 is singular while all proper subsets are
    bent.  It prevents replacing the extension test by a pairwise-only test.
    """
    base = mm_signs(list(range(8)), 0)
    masks = []
    for ui, vj in ((0, 1), (1, 2), (2, 0)):
        masks.append(sum((((x >> ui) & 1) * ((x >> (3 + vj)) & 1)) << x
                         for x in range(POINTS)))
    signs = np.array([[base[x] * (-1 if sum(((z >> j) & 1) * ((masks[j] >> x) & 1)
                                             for j in range(3)) & 1 else 1)
                       for x in range(POINTS)] for z in range(8)], dtype=np.int16)
    cpu_raw = signs @ hmat
    gpu_raw, gpu_seconds = mps_candidates(signs, h_mps)
    if not np.array_equal(gpu_raw, cpu_raw.astype(np.float32)):
        raise AssertionError("pairwise trap MPS/CPU mismatch")
    legal = np.all(np.abs(cpu_raw) == TARGET, axis=1)
    if legal.tolist() != [True] * 7 + [False]:
        raise AssertionError(("pairwise trap unexpected legality", legal.tolist()))
    if set(np.abs(cpu_raw[7]).tolist()) != {0, 16}:
        raise AssertionError(("pairwise trap triple output", cpu_raw[7].tolist()))
    return {"base": "(-1)^(u dot v)", "masks": ["u0v1", "u1v2", "u2v0"],
            "proper_subsets_bent": bool(np.all(legal[:7])), "triple_rejected": not bool(legal[7]),
            "triple_absolute_values": sorted(set(map(int, np.abs(cpu_raw[7])))),
            "mps_cpu_entrywise_equal": True, "mps_float32_transfer_seconds": gpu_seconds,
            "tolerance": 0}


def q_candidates(legal_flat_masks, c_pivots):
    by_q = {}
    for mask in legal_flat_masks:
        q = reduce_by(mask, c_pivots)
        if q:
            by_q.setdefault(q, mask)
    return by_q


def true_walsh_legality_oracle(base):
    """Cached exact test for an arbitrary 64-bit relative sign mask.

    Crucially, a legal XOR of two affine-flat indicators need not itself be an
    affine-flat indicator.  The old debug implementation made that invalid
    restriction; all span extensions must instead call this oracle.
    """
    cache = {0: True}
    def legal(mask):
        answer = cache.get(mask)
        if answer is None:
            signed = [(-value if (mask >> x) & 1 else value) for x, value in enumerate(base)]
            answer = all(abs(value) == TARGET for value in fwht(signed))
            cache[mask] = answer
        return answer
    return legal, cache


def v_slice_positive_baseline(base, hmat, h_mps, c_pivots):
    """Eight flats {v=a}; their full W span is legal for every MM seed."""
    slices = [sum(1 << (u | (a << 3)) for u in range(8)) for a in range(8)]
    signs = np.array([[base[x] * (-1 if sum(((z >> j) & 1) * ((slices[j] >> x) & 1)
                                             for j in range(8)) & 1 else 1)
                       for x in range(POINTS)] for z in range(256)], dtype=np.int16)
    cpu_raw = signs @ hmat
    gpu_raw, gpu_seconds = mps_candidates(signs, h_mps)
    if not np.array_equal(gpu_raw, cpu_raw.astype(np.float32)):
        raise AssertionError("v=a MPS/CPU mismatch")
    if not np.all(np.abs(cpu_raw) == TARGET):
        raise AssertionError("a v-slice sign pattern was not bent")
    integer_verify(base, slices)
    total_rank = rank(list(c_pivots.values()) + slices)
    qrank = rank(list(q_candidates(set(slices), c_pivots)))
    if rank(slices) != 8 or total_rank != 11 or qrank != 4:
        raise AssertionError((rank(slices), total_rank, qrank))
    return {"description": "all eight 3-flats {v=a}, equivalently arbitrary h changes",
            "flat_span_dimension": 8, "dimC_plus_flat_span": 11,
            "extra_dimension_mod_C": 4, "all_W_slices_integer_verified": 256,
            "mps_cpu_entrywise_equal": True, "mps_float32_transfer_seconds": gpu_seconds,
            "tolerance": 0}


def extendable(v, span, legal):
    return all(legal(v ^ w) for w in span)


def random_greedy(candidates, legal, c_pivots, rng, trials, beam_width):
    """Repeated random order and a bounded beam, testing the whole current span."""
    qs = list(candidates)
    best = []
    attempts = 0
    max_examined = 0
    for _ in range(trials):
        rng.shuffle(qs)
        basis, span = [], {0}
        for q in qs:
            v = candidates[q]
            if rank(basis + [q]) == len(basis) + 1 and extendable(v, span, legal):
                basis.append(q)
                span |= {v ^ w for w in list(span)}
        attempts += 1
        max_examined = max(max_examined, len(span))
        if len(basis) > len(best):
            best = basis
    # Beam states are (quotient basis, actual W span); a state is keyed by its span.
    states = [([], {0})]
    for _depth in range(max(0, 9 - len(best))):
        next_states = []
        seen = set()
        for basis, span in states:
            order = qs[:]
            rng.shuffle(order)
            for q in order:
                if rank(basis + [q]) != len(basis) + 1:
                    continue
                v = candidates[q]
                if not extendable(v, span, legal):
                    continue
                new_span = span | {v ^ w for w in span}
                key = tuple(sorted(new_span))
                if key not in seen:
                    seen.add(key)
                    next_states.append((basis + [q], new_span))
        if not next_states:
            break
        # Prefer higher-dimensional spans, then states with most legal extensions.
        rng.shuffle(next_states)
        def score(state):
            basis, span = state
            potential = sum(rank(basis + [q]) == len(basis) + 1 and
                            extendable(candidates[q], span, legal) for q in qs)
            return (len(basis), potential)
        next_states.sort(key=score, reverse=True)
        states = next_states[:beam_width]
        for basis, span in states:
            attempts += 1
            max_examined = max(max_examined, len(span))
            if len(basis) > len(best):
                best = basis[:]
    actual = [candidates[q] for q in best]
    return best, actual, {"random_orders": trials, "states_examined": attempts,
                          "largest_W_slice_count": max_examined, "beam_width": beam_width}


def integer_verify(base, lifts):
    for z in range(1 << len(lifts)):
        signed = base.copy().tolist()
        for j, mask in enumerate(lifts):
            if (z >> j) & 1:
                for x in range(POINTS):
                    if (mask >> x) & 1:
                        signed[x] = -signed[x]
        raw = fwht(signed)
        if any(abs(v) != TARGET for v in raw):
            raise AssertionError((z, raw))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seeds", type=int, default=128)
    parser.add_argument("--orders", type=int, default=96)
    parser.add_argument("--beam-width", type=int, default=48)
    parser.add_argument("--rng-seed", type=int, default=20261117)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    torch.set_num_threads(1)
    if not torch.backends.mps.is_available():
        raise RuntimeError("MPS unavailable; explicit GPU task must not silently fall back")
    hmat = walsh_matrix()
    flat_masks, flat_bits = enumerate_affine_3flats()
    c_pivots = c_space()
    h_mps = torch.from_numpy(hmat.astype(np.float32)).to("mps")
    torch.mps.synchronize()
    trap = pairwise_trap(hmat, h_mps)
    rng = random.Random(args.rng_seed)
    rows, global_best = [], None
    total_cpu = total_gpu = 0.0
    for seed_index in range(args.seeds):
        permutation = list(range(8))
        rng.shuffle(permutation)
        h = rng.randrange(256)
        base = mm_signs(permutation, h)
        base_raw = fwht(base.tolist())
        if any(abs(v) != TARGET for v in base_raw):
            raise AssertionError("MM seed was not bent")
        cpu_began = time.perf_counter()
        signs, cpu_raw, cpu_legal = candidate_masks(base, flat_bits, hmat)
        cpu_elapsed = time.perf_counter() - cpu_began
        gpu_raw, gpu_elapsed = mps_candidates(signs, h_mps)
        if not np.array_equal(gpu_raw, cpu_raw.astype(np.float32)):
            mismatch = int(np.count_nonzero(gpu_raw != cpu_raw.astype(np.float32)))
            raise AssertionError(("MPS/CPU entry mismatch", seed_index, mismatch))
        gpu_legal = np.all(np.abs(gpu_raw) == float(TARGET), axis=1)
        if not np.array_equal(gpu_legal, cpu_legal):
            raise AssertionError(("MPS/CPU legality mismatch", seed_index))
        legal_masks = {flat_masks[i] for i, ok in enumerate(cpu_legal) if ok}
        candidates = q_candidates(legal_masks, c_pivots)
        legal, legality_cache = true_walsh_legality_oracle(base)
        q_basis, actual_basis, search = random_greedy(candidates, legal, c_pivots,
                                                        rng, args.orders, args.beam_width)
        integer_verify(base, actual_basis)
        positive = v_slice_positive_baseline(base, hmat, h_mps, c_pivots)
        row = {"seed_index": seed_index, "permutation": permutation, "h_bits": h,
               "base_inverse_walsh_abs": TARGET, "affine_flats_tested": len(flat_masks),
               "legal_affine_flat_masks": len(legal_masks),
               "nontrivial_quotient_candidates": len(candidates),
               "extra_basis_dimension": len(actual_basis),
               "dimC": N + 1, "dimV": N + 1 + len(actual_basis),
               "affine_support_rank": N + len(actual_basis),
               "quotient_basis_masks": q_basis, "affine_flat_basis_masks": actual_basis,
               "integer_slices_verified": 1 << len(actual_basis),
               "exact_arbitrary_xor_masks_checked": len(legality_cache),
               "v_slice_positive_baseline": positive,
               "cpu_integer_seconds": cpu_elapsed, "mps_float32_transfer_seconds": gpu_elapsed,
               "mps_cpu_entrywise_equal": True, "tolerance": 0, "search": search}
        rows.append(row)
        total_cpu += cpu_elapsed
        total_gpu += gpu_elapsed
        if global_best is None or (row["dimV"], row["nontrivial_quotient_candidates"]) > (global_best["dimV"], global_best["nontrivial_quotient_candidates"]):
            global_best = row
    torch.mps.empty_cache()
    result = {"method": "MPS_float32_Walsh_screen_then_exact_integer_F2_span_verification",
              "scope": "random_Maiorana_McFarland_6_variable_bent_seeds_and_all_3_dimensional_affine_flat_sign_flips",
              "claim_limit": "finite restricted-family search; no general extremal conclusion",
              "points": POINTS, "affine_flats": len(flat_masks), "expected_affine_flats": 11160,
              "seeds": args.seeds, "orders_per_seed": args.orders, "beam_width": args.beam_width,
              "rng_seed": args.rng_seed, "device": "mps:0", "mps_available": True,
              "torch": torch.__version__, "python": sys.version, "platform": platform.platform(),
              "cpu_threads": torch.get_num_threads(), "arithmetic": {"GPU": "float32, exact integer-valued range [-64,64] checked entrywise", "verification": "Python integer FWHT", "tolerance": 0},
              "timing_seconds": {"cpu_integer_total": total_cpu, "mps_transfer_and_compute_total": total_gpu,
                                 "per_seed_cpu_integer_mean": total_cpu / args.seeds,
                                 "per_seed_mps_mean": total_gpu / args.seeds},
              "pairwise_trap": trap, "best": global_best, "rows": rows,
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": "passed", "best_dimV": global_best["dimV"],
                      "best_extra_dimension": global_best["extra_basis_dimension"],
                      "seeds": args.seeds, "mps_total_seconds": total_gpu,
                      "cpu_total_seconds": total_cpu}))


if __name__ == "__main__":
    main()
