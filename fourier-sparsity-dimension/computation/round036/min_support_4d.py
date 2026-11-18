#!/usr/bin/env python3
"""Exact minimum actual supports for common switch codes of the standard four-bit bent phase."""

import argparse
from collections import Counter
from functools import reduce
import hashlib
from itertools import combinations, product
import json
import operator
from pathlib import Path
import platform
import runpy
import sys
import time


def span(basis):
    vectors = {0}
    for v in basis:
        vectors |= {v ^ w for w in vectors}
    return vectors


def basis(space):
    generators = []
    covered = {0}
    for v in sorted(space):
        if v not in covered:
            generators.append(v)
            covered |= {v ^ w for w in covered}
    assert covered == space
    return generators


def contained_subspaces(available, max_dimension):
    """Enumerate by adjoining a complete coset to each lower-dimensional space."""
    levels = [{frozenset({0})}]
    for _ in range(max_dimension):
        extensions = set()
        for space in levels[-1]:
            for v in available - space:
                coset = {v ^ w for w in space}
                if coset <= available:
                    extensions.add(space | coset)
        levels.append(extensions)
    return levels


def enumerate_supports(affine, quadratic, dimensions, legal, f, direct_walsh):
    affine_spaces = contained_subspaces(affine, dimensions[-1])
    quadratic_spaces = contained_subspaces(quadratic, dimensions[-1])
    kernels = []
    for level in affine_spaces:
        records = []
        for space in sorted(level, key=sorted):
            representatives = sorted({min(v ^ w for w in space) for v in affine})
            records.append((basis(space), reduce(operator.or_, space, 0), representatives))
        kernels.append(records)
    projected_bases = [[basis(s) for s in sorted(level, key=sorted)] for level in quadratic_spaces]

    rows = {}
    for dim in dimensions:
        weight_counts = Counter()
        breakdown = {}
        best_weight = len(POINTS) + 1
        best_basis = None
        best_support = None
        for projected_dim in range(dim + 1):
            kernel_dim = dim - projected_dim
            counts = Counter()
            for kernel_basis, kernel_support, cosets in kernels[kernel_dim]:
                for projected_basis in projected_bases[projected_dim]:
                    for shifts in product(cosets, repeat=projected_dim):
                        lifts = [q ^ a for q, a in zip(projected_basis, shifts)]
                        support = reduce(operator.or_, lifts, kernel_support)
                        weight = support.bit_count()
                        counts[weight] += 1
                        if weight < best_weight:
                            best_basis = kernel_basis + lifts
                            best_support = support
                            best_weight = weight
            expected = (
                len(kernels[kernel_dim])
                * len(projected_bases[projected_dim])
                * (len(affine) >> kernel_dim) ** projected_dim
            )
            assert sum(counts.values()) == expected
            breakdown[projected_dim] = dict(counts)
            weight_counts.update(counts)

        code = sorted(span(best_basis))
        assert len(code) == 1 << dim and set(code) <= legal
        assert best_support == reduce(operator.or_, code, 0)
        spectra = {hex(v): direct_walsh(f ^ v) for v in code}
        assert all(abs(value) == 4 for spectrum in spectra.values() for value in spectrum)
        rows[dim] = {
            "subspace_count": sum(weight_counts.values()),
            "support_weight_histogram": dict(weight_counts),
            "support_histogram_by_projection_dim": breakdown,
            "minimum_support_weight": best_weight,
            "minimizing_subspace_count": weight_counts[best_weight],
            "witness": {
                "basis_masks": [hex(v) for v in best_basis],
                "support_mask": hex(best_support),
                "support_points": [v for v in POINTS if (best_support >> v) & 1],
                "all_walsh_spectra_by_switch_mask": spectra,
            },
        }
    return rows, [len(level) for level in affine_spaces], [len(level) for level in quadratic_spaces]


POINTS = range(16)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("round035_source", type=Path, help="R35 exhaust_4d_local.py")
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    started = time.perf_counter()

    r35 = runpy.run_path(str(args.round035_source))
    truth, fwht = r35["truth"], r35["fwht"]
    affine_basis = [(1 << len(POINTS)) - 1] + [
        truth((x >> bit) & 1 for x in POINTS) for bit in range(4)
    ]
    quadratic_basis = [
        truth(((x >> i) & (x >> j)) & 1 for x in POINTS)
        for i, j in combinations(range(4), 2)
    ]
    affine = span(affine_basis)
    homogeneous_quadratic = span(quadratic_basis)
    assert len(affine) == 32 and len(homogeneous_quadratic) == 64
    assert affine & homogeneous_quadratic == {0}
    f = quadratic_basis[0] ^ quadratic_basis[-1]
    assert f == 0x7888

    bent = {mask for mask in range(1 << len(POINTS)) if all(abs(v) == 4 for v in fwht(mask))}
    legal = {f ^ phase for phase in bent}
    quadratic = legal & homogeneous_quadratic
    assert len(bent) == 896 and legal == {q ^ a for q in quadratic for a in affine}
    rows, affine_sizes, quadratic_sizes = enumerate_supports(
        affine, quadratic, (1, 2, 3, 4), legal, f, r35["direct_walsh"]
    )
    result = {
        "scope": "Standard four-variable bent phase f=x0*x1+x2*x3 over F2.",
        "base_phase_truth_mask": hex(f),
        "python": platform.python_version(),
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "round035_source_sha256": hashlib.sha256(args.round035_source.read_bytes()).hexdigest(),
        "arithmetic": "Exact integers, complete enumeration, no RNG or tolerance.",
        "truth_tables_checked": 1 << len(POINTS),
        "legal_directions_including_zero": len(legal),
        "affine_subspace_counts_for_dimensions_0_to_4": affine_sizes,
        "legal_quadratic_subspace_counts_for_dimensions_0_to_4": quadratic_sizes,
        "results_by_dim": rows,
        "elapsed_seconds": time.perf_counter() - started,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({dim: row["minimum_support_weight"] for dim, row in rows.items()}))


if __name__ == "__main__":
    main()
