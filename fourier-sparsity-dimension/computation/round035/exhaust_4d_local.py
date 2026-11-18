#!/usr/bin/env python3
"""Exact four-variable count of bent eight-point blocks and local switch spaces."""

import argparse
from collections import Counter, defaultdict
import hashlib
from itertools import combinations
import json
from pathlib import Path
import platform
import sys
import time


POINTS = range(16)


def truth(values):
    return sum(value << x for x, value in enumerate(values))


def support(mask):
    return [x for x in POINTS if (mask >> x) & 1]


def fwht(mask):
    spectrum = [1 - 2 * ((mask >> x) & 1) for x in POINTS]
    width = 1
    while width < len(POINTS):
        for start in range(0, len(POINTS), 2 * width):
            for left in range(start, start + width):
                right = left + width
                a, b = spectrum[left], spectrum[right]
                spectrum[left] = a + b
                spectrum[right] = a - b
        width *= 2
    return spectrum


def direct_walsh(mask):
    return [
        sum((-1) ** (((mask >> x) & 1) ^ ((x & xi).bit_count() % 2)) for x in POINTS)
        for xi in POINTS
    ]


def gf2_rank(rows):
    pivots = {}
    for row in rows:
        while row:
            pivot = row.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = row
                break
            row ^= pivots[pivot]
    return len(pivots)


def bent_quadratics():
    """Independent calculation of bent truth tables via their polar forms."""
    pairs = list(combinations(range(4), 2))
    result = set()
    for selected in range(1 << len(pairs)):
        edges = [edge for i, edge in enumerate(pairs) if (selected >> i) & 1]
        rows = [0] * 4
        for a, b in edges:
            rows[a] ^= 1 << b
            rows[b] ^= 1 << a
        if gf2_rank(rows) != 4:
            continue
        quadratic = truth(sum(((x >> a) & (x >> b) & 1) for a, b in edges) % 2 for x in POINTS)
        for lin in POINTS:
            affine = truth((lin & x).bit_count() % 2 for x in POINTS)
            result.add(quadratic ^ affine)
            result.add(quadratic ^ affine ^ 0xFFFF)
    return result


def spaces_modulo_block(block, legal):
    """Enumerate all subspaces in the quotient by the block indicator."""
    anchor = block & -block
    available = {
        g for g in legal if g & ~block == 0 and g & anchor == 0 and g ^ block in legal
    }
    families = [{frozenset({0})}]
    while families[-1]:
        next_family = set()
        for space in families[-1]:
            for candidate in available - space:
                coset = {candidate ^ g for g in space}
                if coset <= available:
                    next_family.add(space | coset)
        if not next_family:
            break
        families.append(next_family)
    return families


def witness(block, space, base):
    span = {0}
    basis = [block]
    for g in sorted(space):
        if g not in span:
            basis.append(g)
            span |= {v ^ g for v in span}
    assert span == space
    full_space = sorted(span | {g ^ block for g in span})
    assert len(full_space) == 1 << len(basis)
    assert all(g & ~block == 0 for g in full_space)
    spectra = {hex(g): direct_walsh(base ^ g) for g in full_space}
    assert all(abs(value) == 4 for spectrum in spectra.values() for value in spectrum)
    return {
        "dim": len(basis),
        "block_mask": hex(block),
        "block_support": support(block),
        "affine_span_dimension": gf2_rank(x ^ support(block)[0] for x in support(block)),
        "basis_masks": [hex(g) for g in basis],
        "basis_supports": [support(g) for g in basis],
        "all_walsh_spectra_by_switch_mask": spectra,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    started = time.perf_counter()

    f = truth((((x >> 0) & (x >> 1)) ^ ((x >> 2) & (x >> 3))) & 1 for x in POINTS)
    bent = {mask for mask in range(1 << len(POINTS)) if all(abs(value) == 4 for value in fwht(mask))}
    assert len(bent) == 896 and bent == bent_quadratics()
    legal = {f ^ phase for phase in bent}
    blocks = sorted(g for g in legal if g.bit_count() == 8)
    space_counts = defaultdict(Counter)
    block_counts = defaultdict(Counter)
    witnesses = {}
    for block in blocks:
        block_points = support(block)
        affine_rank = gf2_rank(x ^ block_points[0] for x in block_points)
        assert affine_rank in (3, 4)
        kind = "affine" if affine_rank == 3 else "nonaffine"
        families = spaces_modulo_block(block, legal)
        block_counts[kind][len(families)] += 1
        for dim, family in enumerate(families, start=1):
            space_counts[kind][dim] += len(family)
        if kind not in witnesses or witnesses[kind]["dim"] < len(families):
            canon_space = min(families[-1], key=lambda space: sorted(space))
            witnesses[kind] = witness(block, canon_space, f)

    result = {
        "base_phase_formula": "x0*x1+x2*x3 over F2, least significant bit x0",
        "base_phase_truth_mask": hex(f),
        "python": platform.python_version(),
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "arithmetic": "Exact integer arithmetic, exhaustive search, no RNG or tolerance.",
        "truth_tables_checked": 1 << len(POINTS),
        "bent_truth_tables": len(bent),
        "independent_bent_quadratic_check": True,
        "eight_point_legal_blocks": len(blocks),
        "blocks_by_kind_and_max_local_dim": block_counts,
        "block_space_pairs_by_kind_and_dim": space_counts,
        "max_dim_witnesses": witnesses,
        "elapsed_seconds": time.perf_counter() - started,
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"blocks": block_counts, "spaces": space_counts}, sort_keys=True))


if __name__ == "__main__":
    main()
