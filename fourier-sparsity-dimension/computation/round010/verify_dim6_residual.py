#!/usr/bin/env python3
"""Exact audit of the 13-dimensional global-plus-residual D6 construction.

V = RM(1,6) + {1_R ell(x,b): ell in RM(1,5)}.
This checks all 8192 F+v for bentness and all affine 3-flats for a common
minimum-word flip: F|A and every residual generator |A must be affine.
"""
import argparse, hashlib, json, sys, time
from pathlib import Path


def fwht(a):
    a = a[:]
    h = 1
    while h < len(a):
        for p in range(0, len(a), 2 * h):
            for i in range(p, p + h):
                left, right = a[i], a[i + h]
                a[i], a[i + h] = left + right, left - right
        h *= 2
    return a


def f0(x):
    return ((x & 1) * ((x >> 3) & 1) ^ ((x >> 1) & 1) * ((x >> 2) & 1)) & 1


def vfun(x):
    return ((x & 1) * ((x >> 2) & 1)
            ^ ((x >> 1) & 1) * ((x >> 2) & 1)
            ^ ((x >> 1) & 1) * ((x >> 3) & 1)) & 1


def Fval(p):
    x, b, z = p & 15, (p >> 4) & 1, (p >> 5) & 1
    return f0(x) ^ (z & (b ^ vfun(x)))


def affine_values(values):
    """Whether a truth table on F_2^3 (in standard binary order) is affine."""
    a = list(values)
    # Boolean ANF Möbius transform.
    for bit in range(3):
        for x in range(8):
            if (x >> bit) & 1:
                a[x] ^= a[x ^ (1 << bit)]
    return all(a[x] == 0 for x in range(8) if x.bit_count() >= 2)


def rank3(vs):
    pivots = []
    for x in vs:
        y = x
        for p in pivots:
            y = min(y, y ^ p)
        if not y:
            return False
        pivots.append(y)
        pivots.sort(reverse=True)
        # Maintain a simple reduced independent list.
        for i in range(len(pivots)):
            for j in range(i + 1, len(pivots)):
                if (pivots[i] ^ pivots[j]) < pivots[i]:
                    pivots[i] ^= pivots[j]
        pivots.sort(reverse=True)
    return True


def subspace(vs):
    return frozenset(
        ((combo >> 0 & 1) * vs[0]) ^ ((combo >> 1 & 1) * vs[1]) ^ ((combo >> 2 & 1) * vs[2])
        for combo in range(8)
    )


def linear_3spaces():
    # Store one ordered independent basis for each of the 1395 subspaces.
    spaces = {}
    vectors = range(1, 64)
    for a in vectors:
        for b in range(a + 1, 64):
            for c in range(b + 1, 64):
                basis = (a, b, c)
                if rank3(basis):
                    S = subspace(basis)
                    spaces.setdefault(S, basis)
    return spaces


def coset_reps(S):
    seen, reps = set(), []
    for p in range(64):
        C = frozenset(p ^ s for s in S)
        if C not in seen:
            seen.add(C)
            reps.append(p)
    return reps


def values_on_flat(fun, origin, basis):
    return [fun(origin ^ ((w & 1) * basis[0]) ^ (((w >> 1) & 1) * basis[1])
                ^ (((w >> 2) & 1) * basis[2])) for w in range(8)]


def main():
    ap=argparse.ArgumentParser();ap.add_argument("--out",type=Path,required=True);args=ap.parse_args();began=time.monotonic()
    N, q = 64, 8
    F = [Fval(p) for p in range(N)]
    # Basis for C_global=RM(1,6): constant plus six coordinate bits.
    global_basis = [[1] * N] + [[(p >> i) & 1 for p in range(N)] for i in range(6)]
    # Basis of G=1_R RM(1,5)(x,b): constant, x_0,...,x_3,b.
    residual_basis = []
    for coord in range(6):
        arr = []
        for p in range(N):
            x, b = p & 15, (p >> 4) & 1
            in_R = b == vfun(x)
            ell = 1 if coord == 0 else ((x >> (coord - 1)) & 1 if coord <= 4 else b)
            arr.append(int(in_R and ell))
        residual_basis.append(arr)
    basis = global_basis + residual_basis

    # GF(2) rank of function vectors.
    rows = [sum(bit << p for p, bit in enumerate(row)) for row in basis]
    piv = []
    for row in rows:
        y = row
        for p in piv:
            y = min(y, y ^ p)
        if y:
            piv.append(y)
            piv.sort(reverse=True)

    all_bent, amplitudes, weights = True, set(), set()
    for mask in range(1 << len(basis)):
        phase = [F[p] for p in range(N)]
        wordwt = 0
        for j, row in enumerate(basis):
            if (mask >> j) & 1:
                for p in range(N):
                    phase[p] ^= row[p]
        W = fwht([1 if bit == 0 else -1 for bit in phase])
        aset = {abs(y) for y in W}
        all_bent &= (aset == {q})
        amplitudes |= aset
        # Recompute codeword weight from its value, avoiding phase/base mixing.
        word = [0] * N
        for j, row in enumerate(basis):
            if (mask >> j) & 1:
                word = [x ^ y for x, y in zip(word, row)]
        weights.add(sum(word))

    spaces = linear_3spaces()
    valid = []
    all_flat_count = 0
    for S, bas in spaces.items():
        for origin in coset_reps(S):
            all_flat_count += 1
            # All global-affine basis restrictions are automatically affine.
            tests = [Fval] + [
                (lambda p, row=row: row[p]) for row in residual_basis
            ]
            if all(affine_values(values_on_flat(fun, origin, bas)) for fun in tests):
                flat = sorted(origin ^ s for s in S)
                valid.append({
                    "origin": origin, "basis": list(bas), "points": flat,
                    "F_restriction": values_on_flat(Fval, origin, bas),
                    "residual_restrictions": [
                        values_on_flat(lambda p, row=row: row[p], origin, bas)
                        for row in residual_basis
                    ],
                })

    residual_weight_histogram={}; partial_supports=set(); residual_masks=[]
    for mask in range(64):
        word=[0]*N
        for j,row in enumerate(residual_basis):
            if mask>>j&1: word=[a^b for a,b in zip(word,row)]
        residual_masks.append(sum(bit<<p for p,bit in enumerate(word)))
        wt=sum(word);residual_weight_histogram[wt]=residual_weight_histogram.get(wt,0)+1
        partial=[(1-2*(F[p]^word[p])) if ((p>>4)&1)==vfun(p&15) else 0 for p in range(N)]
        W=fwht(partial);support=tuple(i for i,x in enumerate(W) if x)
        assert set(abs(x) for x in W if x)=={8}
        partial_supports.add(support)
    assert len(set(residual_masks))==64 and len(partial_supports)==1 and next(iter(partial_supports))==tuple(range(32))
    assert 8 not in residual_weight_histogram
    assert len(piv)==13 and all_bent and 8 not in weights and len(spaces)==1395 and all_flat_count==11160
    assert len(valid)==40
    out = {
        "status":"passed","arithmetic":"integer/F2; zero tolerance",
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "python":sys.version,"elapsed_seconds":time.monotonic()-began,
        "F_truthmask_hex":hex(sum(bit<<p for p,bit in enumerate(F))),
        "global_basis_masks_hex":[hex(sum(bit<<p for p,bit in enumerate(row))) for row in global_basis],
        "residual_basis_masks_hex":[hex(sum(bit<<p for p,bit in enumerate(row))) for row in residual_basis],
        "residual_weight_histogram":residual_weight_histogram,
        "residual_common_support":list(next(iter(partial_supports))),
        "q": q,
        "basis_dimension_expected": 13,
        "basis_rank": len(piv),
        "candidate_space_size": 1 << len(piv),
        "all_8192_F_plus_v_bent": all_bent,
        "Walsh_absolute_value_set_across_slices": sorted(amplitudes),
        "V_weight_set": sorted(weights),
        "V_has_q_weight": q in weights,
        "linear_3space_count": len(spaces),
        "affine_3flat_count": all_flat_count,
        "common_minimum_flip_flat_count": len(valid),
        "common_minimum_flip_flats": valid,
    }
    with args.out.open("w") as f:
        json.dump(out, f, indent=2, sort_keys=True)
    print(json.dumps({key:value for key,value in out.items() if key not in ("common_minimum_flip_flats","residual_common_support")}))


if __name__ == "__main__":
    main()
