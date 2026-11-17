#!/usr/bin/env python3
"""Exact common-MM-splitting audit for round005 dimV=15 certificates.

For every 3-space U <= F_2^6, test whether the seed phase q and every one
of the eight stored relative-sign masks are affine on every U-coset.  The
criterion is D_a D_b f = 0 for the three pairs in a basis (a,b,c) of U and
all 64 base points.  It is exact over F_2.
"""
import argparse
import hashlib
import itertools
import json
import platform
import sys
import time
from pathlib import Path

POINTS, DIM = 64, 6

def xor_span_mask(subspace_mask, direction):
    members = [x for x in range(POINTS) if subspace_mask >> x & 1]
    out = subspace_mask
    for x in members:
        out |= 1 << (x ^ direction)
    return out

def rref_basis(subspace_mask):
    piv = [0] * DIM
    for x in range(1, POINTS):
        if not (subspace_mask >> x & 1):
            continue
        y = x
        for bit in range(DIM - 1, -1, -1):
            if y >> bit & 1:
                if piv[bit]:
                    y ^= piv[bit]
                else:
                    piv[bit] = y
                    for other in range(DIM):
                        if other != bit and (piv[other] >> bit & 1):
                            piv[other] ^= y
                    break
    out = tuple(x for x in piv if x)
    if len(out) != 3:
        raise AssertionError((subspace_mask, out))
    return out

def all_three_subspaces():
    levels = [{1}]  # bit x records membership of x; {0} is mask 1.
    for dim in range(3):
        following = set()
        for space in levels[-1]:
            for direction in range(1, POINTS):
                if not (space >> direction & 1):
                    following.add(xor_span_mask(space, direction))
        levels.append(following)
    expected = (1, 63, 651, 1395)
    if tuple(map(len, levels)) != expected:
        raise AssertionError(tuple(map(len, levels)))
    return [(rref_basis(space), space) for space in sorted(levels[3])]

def seed_phase(permutation, h_bits):
    return [(((x & 7) & permutation[x >> 3]).bit_count() ^ ((h_bits >> (x >> 3)) & 1)) & 1
            for x in range(POINTS)]

def mask_phase(mask):
    return [(mask >> x) & 1 for x in range(POINTS)]

def affine_on_all_cosets_by_differences(values, basis):
    for a, b in itertools.combinations(basis, 2):
        for x in range(POINTS):
            if values[x] ^ values[x ^ a] ^ values[x ^ b] ^ values[x ^ a ^ b]:
                return False
    return True

def load_certificates(paths):
    certificates = []
    for path in paths:
        raw = path.read_bytes()
        data = json.loads(raw)
        for row_index, row in enumerate(data["rows"]):
            if row.get("dimV") != 15 or len(row["affine_flat_basis_masks"]) != 8:
                raise AssertionError((path, row_index, row.get("dimV")))
            certificates.append({"source_file": str(path), "source_sha256": hashlib.sha256(raw).hexdigest(),
                                 "row_index": row_index, "seed_index": row["seed_index"],
                                 "permutation": row["permutation"], "h_bits": row["h_bits"],
                                 "relative_generators": row["affine_flat_basis_masks"]})
    return certificates

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, action="append", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    spaces = all_three_subspaces()
    certs = load_certificates(args.input)
    rows = []
    for cert in certs:
        functions = [seed_phase(cert["permutation"], cert["h_bits"])]
        functions.extend(mask_phase(mask) for mask in cert["relative_generators"])
        feasible = []
        for basis, _ in spaces:
            if all(affine_on_all_cosets_by_differences(f, basis) for f in functions):
                feasible.append(list(basis))
        rows.append({**cert, "functions_tested": 9, "valid_U_count": len(feasible),
                     "valid_U_rref_bases": feasible,
                     "first_valid_U_rref_basis": feasible[0] if feasible else None})
    report = {"method": "all_1395_dim3_subspaces_F2_6_second_directional_differences",
              "criterion": "for q and each of 8 relative generators, D_aD_b f=0 for all x and all 3 basis pairs",
              "subspaces_tested": len(spaces), "certificate_count": len(rows), "rows": rows,
              "elapsed_seconds": time.monotonic() - started, "python": sys.version,
              "platform": platform.platform(),
              "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    args.out.write_text(text)
    print(json.dumps({"status": "passed", "certificate_count": len(rows), "subspaces_tested": len(spaces),
                      "valid_U_counts": [r["valid_U_count"] for r in rows], "elapsed_seconds": report["elapsed_seconds"]}))

if __name__ == "__main__":
    main()
