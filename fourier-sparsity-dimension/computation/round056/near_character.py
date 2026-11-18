#!/usr/bin/env python3
"""Detect and reconstruct a signed character flipped on one affine coset."""

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from functools import cache
import hashlib
import json
from pathlib import Path
import platform
import time


def parity(word):
    return word.bit_count() & 1


def points(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def basis(words):
    pivots, result = {}, []
    for word in words:
        x = word
        while x:
            k = x.bit_length() - 1
            if k not in pivots:
                pivots[k] = x
                result.append(word)
                break
            x ^= pivots[k]
    return result


def fwht(values):
    out = list(values)
    step = 1
    while step < len(out):
        for start in range(0, len(out), 2 * step):
            for i in range(start, start + step):
                j = i + step
                a, b = out[i], out[j]
                out[i], out[j] = a + b, a - b
        step *= 2
    return out


def linear_solve(equations):
    echelon = {}
    for q, rhs in equations:
        while q:
            k = q.bit_length() - 1
            if k not in echelon:
                echelon[k] = (q, rhs)
                break
            row, r = echelon[k]
            q ^= row
            rhs ^= r
        if not q and rhs:
            return None
    x = 0
    for k, (row, rhs) in sorted(echelon.items()):
        x |= (rhs ^ parity(x & row)) << k
    return x


def span(rows):
    output = [0]
    for q in rows:
        output.extend(x ^ q for x in output[:])
    return output


def signed_character_value(f, x):
    twist = parity(f["carrier"] & x) ^ f["carrier_negative"]
    normals = f["normal_basis"]
    flipped = normals and all(not parity(row & (x ^ f["affine_origin"])) for row in normals)
    return 1 - 2 * (twist ^ int(bool(flipped)))


def reconstruct(spectrum, ambient_bits, selected_carrier=None):
    """Return a coset representation from exact integer Walsh data, if structurally valid."""
    M = 1 << ambient_bits
    if not spectrum:
        return None
    peak = max(map(abs, spectrum.values()))
    carrier = min(q for q, val in spectrum.items() if abs(val) == peak) if selected_carrier is None else selected_carrier
    if abs(spectrum[carrier]) != peak:
        return None
    neg = int(spectrum[carrier] < 0)
    if len(spectrum) == 1:
        if peak != M:
            return None
        return {"carrier": carrier, "carrier_negative": neg, "normal_basis": [], "affine_origin": 0,
                "affine_spectrum_rank": 0, "linear_spectrum_rank": int(bool(carrier))}
    normal_points = [a ^ carrier for a in sorted(spectrum)]
    rows = basis(normal_points)
    r = len(rows)
    if r < 2 or len(spectrum) != 1 << r:
        return None
    absolute_offpeak = 1 << (ambient_bits - r + 1)
    if peak != M - absolute_offpeak:
        return None
    sign = 1 - 2 * neg
    equation = []
    for q, coeff in spectrum.items():
        if q == carrier:
            continue
        if abs(coeff) != absolute_offpeak:
            return None
        equation.append((q ^ carrier, int(sign * coeff > 0)))
    x0 = linear_solve(equation)
    if x0 is None:
        return None
    assert all(parity(q & x0) == rhs for q, rhs in equation)
    d = len(basis(spectrum))
    assert d == len(basis(rows + [carrier]))
    return {
        "carrier": carrier, "carrier_negative": neg, "normal_basis": rows, "affine_origin": x0,
        "affine_spectrum_rank": r, "linear_spectrum_rank": d,
    }


def normalized_representation_spectrum(record, n):
    """Generate the full integer spectrum of the reconstructed Boolean formula."""
    a = record["carrier"]
    sign = 1 - 2 * record["carrier_negative"]
    rows = record["normal_basis"]
    M = 1 << n
    if not rows:
        return {a: sign * M}
    step = 1 << (n - len(rows) + 1)
    W = {a: sign * (M - step)}
    for q in span(rows)[1:]:
        W[q ^ a] = -sign * step * (1 - 2 * parity(q & record["affine_origin"]))
    return W


@cache
def geom_catalog():
    # This construction never computes a Walsh transform.
    total = 1 << 16
    chars = [sum(parity(a & x) << x for x in range(16)) for a in range(16)]
    affines = []
    for mask in range(1, total):
        size = mask.bit_count()
        if size & (size - 1):
            continue
        a, *tail = points(mask)
        dim = len(basis(x ^ a for x in tail))
        if size == 1 << dim:
            affines.append((dim, mask))
    assert Counter(dim for dim, _ in affines) == {0: 16, 1: 120, 2: 140, 3: 30, 4: 1}
    words = Counter({char ^ sign: 1 for char in chars for sign in (0, total - 1)})
    for dim, coset in affines:
        if dim >= 3:  # Hyperplane and whole-space flips already give signed characters.
            assert all((char ^ sign ^ coset) in words for char in chars for sign in (0, total - 1))
            continue
        for char in chars:
            for sign in (0, total - 1):
                words[char ^ sign ^ coset] += 1
    assert len(words) == 5504
    return words, chars


def audit_truth_range(interval):
    catalog, chars = geom_catalog()
    lo, hi = interval
    records, rank_counts, norm_counts = {}, Counter(), Counter()
    peak_counts = Counter()
    near, boundary, max_carrier_checks = 0, 0, 0
    for mask in range(lo, hi):
        spectrum = {i: v for i, v in enumerate(fwht(1 - 2 * ((mask >> x) & 1) for x in range(16))) if v}
        numerator = sum(map(abs, spectrum.values()))
        rep = reconstruct(spectrum, 4)
        assert (rep is not None) == (mask in catalog) == (numerator < 48)
        boundary += numerator == 48
        norm_counts[numerator] += 1
        if rep is None:
            continue
        near += 1
        assert normalized_representation_spectrum(rep, 4) == spectrum
        r = rep["affine_spectrum_rank"]
        d = rep["linear_spectrum_rank"]
        largest = max(map(abs, spectrum.values()))
        tied = [a for a, coeff in spectrum.items() if abs(coeff) == largest]
        assert len(tied) == catalog[mask]
        for carrier in tied:
            specialized = reconstruct(spectrum, 4, carrier)
            truth = sum(int(signed_character_value(specialized, x) < 0) << x for x in range(16))
            assert truth == mask
            if r:
                flip = mask ^ chars[carrier] ^ (65535 if specialized["carrier_negative"] else 0)
                assert flip.bit_count() == (1 << (4 - r))
            max_carrier_checks += 1
        rank_counts[r, d, numerator] += 1
        peak_counts[r, len(tied)] += 1
        records.setdefault((r, d), {"truth_mask": hex(mask), "norm": str(Fraction(numerator, 16)), **rep})
    return near, boundary, max_carrier_checks, records, rank_counts, norm_counts, peak_counts


def compressed_lifts():
    results = []
    for n, r in [(8, 2), (9, 3), (12, 4), (15, 5), (40, 6), (64, 8), (70, 10)]:
        normals = [(1 << i) | (1 << ((i * 3 + r + 1) % (n - r) + r)) for i in range(r)]
        assert len(basis(normals)) == r
        outside = 1 << (n - 1)
        assert len(basis(normals + [outside])) == r + 1
        for case, carrier in (("in_direction", normals[0]), ("out_of_direction", normals[0] ^ outside)):
            constant = {"carrier": carrier, "carrier_negative": 1, "normal_basis": normals,
                        "affine_origin": ((1 << (n - 1)) | 5)}
            W = normalized_representation_spectrum(constant, n)
            rep = reconstruct(W, n)
            assert normalized_representation_spectrum(rep, n) == W
            assert rep["linear_spectrum_rank"] == r + int(case == "out_of_direction")
            # Testing one point from every joint character/coset atom covers the whole cube.
            Q = basis(normals + [carrier])
            d = len(Q)
            for para in range(1 << d):
                x = linear_solve([(row, (para >> i) & 1) for i, row in enumerate(Q)])
                assert x is not None
                assert signed_character_value(rep, x) == signed_character_value(constant, x)
            results.append({
                "ambient_dimension": n, "offset_case": case, "affine_spectrum_rank": r,
                "linear_spectrum_rank": d, "joint_coordinate_atoms_covered": 1 << d,
                "norm": str(Fraction(sum(map(abs, W.values())), 1 << n)),
                "recovers_every_coefficient": True, "recovered_representation": rep,
            })
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--workers", type=int, default=12)
    args = parser.parse_args()
    start = time.perf_counter()
    counts, norms, peaks = Counter(), Counter(), Counter()
    near = boundary = max_checks = 0
    records = {}
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        sh = [(lo, lo + 4096) for lo in range(0, 65536, 4096)]
        for k, b, t, reps, r, nl, pk in pool.map(audit_truth_range, sh):
            near += k
            boundary += b
            max_checks += t
            counts.update(r)
            norms.update(nl)
            peaks.update(pk)
            for key, record in reps.items():
                records.setdefault(key, record)
    assert near == 5504 and sum(norms.values()) == 65536
    result = {
        "python": platform.python_version(), "workers": args.workers,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "truth_tables": 65536, "strictly_below_3": near, "exactly_3": boundary,
        "geometric_catalog_exactly_same": True,
        "all_maximal_carrier_reconstructions_checked": max_checks,
        "norm_numerators_over_16": dict(sorted(norms.items())),
        "below_3_by_affine_rank_linear_rank_norm": [[a, d, str(Fraction(n, 16)), v]
                                                  for (a, d, n), v in sorted(counts.items())],
        "carrier_multiplicity_by_affine_rank": [[a, k, v] for (a, k), v in sorted(peaks.items())],
        "canonical_finite_witnesses": list(records.values()),
        "higher_environment_exact_lifts": compressed_lifts(),
        "elapsed_seconds": time.perf_counter() - start,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
