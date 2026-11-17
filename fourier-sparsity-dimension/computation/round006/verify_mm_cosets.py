#!/usr/bin/env python3
"""Independent 8-point-coset affine-table verifier for audit_common_mm_split."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path

POINTS = 64

def popcount(x):
    return bin(x).count("1")

def phase_q(permutation, h_bits):
    return [(popcount((x & 7) & permutation[x >> 3]) ^ ((h_bits >> (x >> 3)) & 1)) & 1
            for x in range(POINTS)]

def coords_and_cosets(basis):
    vectors = []
    for bits in range(8):
        v = 0
        for j in range(3):
            if bits >> j & 1: v ^= basis[j]
        vectors.append(v)
    seen, cosets = set(), []
    for r in range(POINTS):
        orbit = tuple(sorted(r ^ u for u in vectors))
        if orbit not in seen:
            seen.add(orbit); cosets.append((r, vectors))
    if len(cosets) != 8: raise AssertionError((basis, len(cosets)))
    return cosets

def span_mask(basis):
    mask = 0
    for bits in range(8):
        x = 0
        for j in range(3):
            if bits >> j & 1: x ^= basis[j]
        mask |= 1 << x
    return mask

def all_three_spaces_by_independent_triples():
    """Separate enumeration: independent triples, deduplicated by their span."""
    spaces = {}
    for a in range(1, 64):
        for b in range(a + 1, 64):
            for c in range(b + 1, 64):
                mask = span_mask((a, b, c))
                if popcount(mask) == 8:
                    spaces.setdefault(mask, (a, b, c))
    if len(spaces) != 1395: raise AssertionError(len(spaces))
    return spaces

def affine_from_tables(values, basis):
    for rep, vectors in coords_and_cosets(basis):
        c = values[rep]
        slopes = [values[rep] ^ values[rep ^ basis[j]] for j in range(3)]
        for bits, vector in enumerate(vectors):
            expected = c
            for j in range(3):
                if bits >> j & 1: expected ^= slopes[j]
            if values[rep ^ vector] != expected:
                return False
    return True

parser = argparse.ArgumentParser()
parser.add_argument("audit", type=Path)
parser.add_argument("--out", type=Path, required=True)
args = parser.parse_args()
raw = args.audit.read_bytes(); data = json.loads(raw)
if data["subspaces_tested"] != 1395 or data["certificate_count"] != 20: raise AssertionError("scope")
checked = []
spaces = all_three_spaces_by_independent_triples()
for row in data["rows"]:
    functions = [phase_q(row["permutation"], row["h_bits"])]
    functions += [[(m >> x) & 1 for x in range(POINTS)] for m in row["relative_generators"]]
    reported = {span_mask(b) for b in row["valid_U_rref_bases"]}
    independently_feasible = {mask for mask, basis in spaces.items()
                              if all(affine_from_tables(f, basis) for f in functions)}
    if independently_feasible != reported:
        raise AssertionError((row["source_file"], row["row_index"]))
    checked.append({"source_file": row["source_file"], "row_index": row["row_index"],
                    "valid_U_count": len(reported)})
args.out.write_text(json.dumps({"status": "passed", "audit_sha256": hashlib.sha256(raw).hexdigest(),
                                "algorithm": "independent_triple_enumeration_plus_explicit_8_point_coset_affine_tables", "rows": checked}, indent=2)+"\n")
