#!/usr/bin/env python3
"""Independent integer acceptance check for a recorded corrected GPU search."""
import argparse
import hashlib
import json
from pathlib import Path

N, POINTS, TARGET = 6, 64, 8


def rank(vectors):
    pivots = {}
    for v in vectors:
        while v:
            p = v.bit_length() - 1
            if p in pivots:
                v ^= pivots[p]
            else:
                pivots[p] = v
                break
    return len(pivots)


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


def base_signs(permutation, h):
    return [1 if not (((x & 7) & permutation[x >> 3]).bit_count() & 1) ^ ((h >> (x >> 3)) & 1)
            else -1 for x in range(POINTS)]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("result", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    raw = args.result.read_bytes()
    data = json.loads(raw)
    c = [(1 << POINTS) - 1] + [sum(((x >> j) & 1) << x for x in range(POINTS)) for j in range(N)]
    if rank(c) != 7:
        raise AssertionError("C rank")
    rows = []
    for row in data["rows"]:
        basis, qbasis = row["affine_flat_basis_masks"], row["quotient_basis_masks"]
        if len(basis) != row["extra_basis_dimension"] or len(qbasis) != len(basis):
            raise AssertionError("stored basis length")
        if rank(qbasis) != len(qbasis) or rank(c + basis) != 7 + len(basis):
            raise AssertionError(("not independent modulo C", row["seed_index"]))
        base = base_signs(row["permutation"], row["h_bits"])
        for z in range(1 << len(basis)):
            signed = [(-base[x] if sum(((z >> j) & 1) * ((basis[j] >> x) & 1)
                                       for j in range(len(basis))) & 1 else base[x])
                      for x in range(POINTS)]
            if any(abs(v) != TARGET for v in fwht(signed)):
                raise AssertionError(("non-bent W slice", row["seed_index"], z))
        rows.append({"seed_index": row["seed_index"], "extra_dimension": len(basis),
                     "all_W_slices": 1 << len(basis)})
    result = {"status": "passed", "arithmetic": "Python integers", "tolerance": 0,
              "input_sha256": hashlib.sha256(raw).hexdigest(), "rows": rows,
              "total_W_slices": sum(x["all_W_slices"] for x in rows),
              "all_bases_independent_mod_C": True,
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": "passed", "rows": len(rows),
                      "total_W_slices": result["total_W_slices"]}))


if __name__ == "__main__":
    main()
