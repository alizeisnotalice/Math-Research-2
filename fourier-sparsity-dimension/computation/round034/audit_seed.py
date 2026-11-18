#!/usr/bin/env python3
"""Reproduce the finite eight-variable R33 seed calculations used in R34."""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import platform
import sys
import time


POINTS = range(256)
PLANES = (
    (0, 2, 4, 6),
    (65, 67, 68, 70),
    (74, 75, 88, 89),
    (104, 105, 120, 121),
)
EXPECTED_SIGNS = (1, -1, -1, 1, -1, 1, 1, -1, 1, -1, -1, 1, -1, 1, 1, -1)


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def parity(value):
    return value.bit_count() % 2


def truth(values):
    return sum(value << point for point, value in enumerate(values))


def signs(mask):
    return [1 - 2 * ((mask >> point) & 1) for point in POINTS]


def project_bits(point, coordinates):
    return sum(((point >> bit) & 1) << j for j, bit in enumerate(coordinates))


def seed_switch(point):
    x0, x1, x2, x3, x4, x5, x6, x7 = [(point >> bit) & 1 for bit in range(8)]
    return (1 ^ x7) & ((x0 & x1 & x3) ^ (x0 & x4 & x5) ^ (x1 & x2 & x6))


def fwht(values):
    spectrum = list(values)
    width = 1
    while width < len(spectrum):
        for start in range(0, len(spectrum), 2 * width):
            for left in range(start, start + width):
                right = left + width
                a, b = spectrum[left], spectrum[right]
                spectrum[left] = a + b
                spectrum[right] = a - b
        width *= 2
    return spectrum


def direct_walsh(values):
    return [
        sum(value * (1 - 2 * parity(x & xi)) for x, value in enumerate(values))
        for xi in range(len(values))
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


def periods(points, universe_size):
    return [t for t in range(universe_size) if {point ^ t for point in points} == points]


def load_seed(data):
    spec = json.loads(data)["solution"]
    data_bits = spec["data_coords"]
    selector_bits = spec["selector_coords"]
    pi, rho = spec["pi"], spec["rho"]
    assert sorted(data_bits + selector_bits) == list(range(8))
    assert sorted(pi) == list(range(16)) and len(rho) == 16 and set(rho) <= {0, 1}
    b = int(spec["f_hex"], 16)
    b_plus_h = int(spec["f_plus_h_hex"], 16)
    h = b ^ b_plus_h
    for x in POINTS:
        z, w = project_bits(x, data_bits), project_bits(x, selector_bits)
        assert ((b >> x) & 1) == parity(z & pi[w]) ^ rho[w]
        assert ((h >> x) & 1) == seed_switch(x)
    return b, b_plus_h, h


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("seed", type=Path, help="R33 mm_completion_result.json")
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    started = time.perf_counter()

    seed_bytes = args.seed.read_bytes()
    b, b_plus_h, h = load_seed(seed_bytes)
    space_size = len(POINTS)
    q = 16
    R = {x for x in POINTS if (h >> x) & 1}
    assert len(R) == 2 * q
    assert {x ^ y for x in R for y in R} == set(range(128))
    assert periods(R, space_size) == [0]

    base_phase = signs(b)
    switched_phase = signs(b_plus_h)
    base_restriction = [value if x in R else 0 for x, value in enumerate(base_phase)]
    base_full = direct_walsh(base_phase)
    switched_full = direct_walsh(switched_phase)
    base_partial = direct_walsh(base_restriction)
    assert all(abs(value) == q for value in base_full + switched_full)
    assert {i for i in POINTS if base_full[i] != switched_full[i]} == {
        i for i in POINTS if base_partial[i]
    }
    assert switched_full == [total - 2 * part for total, part in zip(base_full, base_partial)]

    D = {i for i in POINTS if base_partial[i]}
    T = {i % 128 for i in D}
    assert len(T) == 16 and len(D) == 32 and D == {i ^ t for i in T for t in (0, 128)}
    assert [base_partial[i] for i in sorted(T)] == [q * s for s in EXPECTED_SIGNS]
    assert periods(T, 128) == [0]
    assert periods(D, space_size) == [0, 128]

    plane_rows = []
    common_direction = set(range(128))
    for plane in PLANES:
        direction = {point ^ plane[0] for point in plane}
        assert len(direction) == 4 and {a ^ b for a in direction for b in direction} == direction
        assert gf2_rank(direction) == 2
        common_direction &= direction
        plane_rows.append({"points": list(plane), "direction": sorted(direction)})
    assert Counter(point for plane in PLANES for point in plane) == Counter(T)
    assert common_direction == {0}

    affine_characters = [truth(parity(x & lam) for x in POINTS) for lam in POINTS]
    constant = (1 << space_size) - 1
    basis = [constant] + [affine_characters[1 << i] for i in range(8)] + [h]
    assert gf2_rank(basis) == 10
    restriction_rank = gf2_rank(row & h for row in basis)
    assert restriction_rank == 8

    supported_parameters = []
    phase_count = 0
    for lam, linear in enumerate(affine_characters):
        for c in (0, 1):
            for e in (0, 1):
                v = linear ^ (c * constant) ^ (e * h)
                phase = signs(b ^ v)
                total = fwht(phase)
                part = fwht([value if x in R else 0 for x, value in enumerate(phase)])
                old = (base_full, switched_full)[e]
                assert total == [(1 - 2 * c) * old[xi ^ lam] for xi in POINTS]
                assert all(abs(value) == q for value in total)
                assert part == [(1 - 2 * (c ^ e)) * base_partial[xi ^ lam] for xi in POINTS]
                assert {i for i, value in enumerate(part) if value} == {i ^ lam for i in D}
                assert {abs(value) for value in part} == {0, q}
                assert (v & h) == ((linear ^ ((c ^ e) * constant)) & h)
                if v & ~h == 0:
                    supported_parameters.append({"c": c, "lambda": lam, "e": e})
                phase_count += 1
    assert phase_count == 1024
    assert supported_parameters == [{"c": 0, "lambda": 0, "e": e} for e in (0, 1)]

    result = {
        "scope": "Exact audit of the eight-variable R33 seed and the displayed common code.",
        "input_basename": args.seed.name,
        "input_sha256": sha256(seed_bytes),
        "script_sha256": sha256(Path(__file__).read_bytes()),
        "python": platform.python_version(),
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "arithmetic": "Exact integer arithmetic; complete enumeration; no RNG or tolerance.",
        "b_hex": hex(b),
        "b_plus_h_hex": hex(b_plus_h),
        "h_hex": hex(h),
        "R_points": sorted(R),
        "R_periods": periods(R, space_size),
        "R_xor_set_dimension": 7,
        "partial_support_T": sorted(T),
        "partial_spectrum_T_values": [base_partial[xi] for xi in sorted(T)],
        "partial_support_D": sorted(D),
        "T_periods": periods(T, 128),
        "D_periods": periods(D, space_size),
        "T_plane_partition": plane_rows,
        "plane_common_direction": sorted(common_direction),
        "code_dimension": gf2_rank(basis),
        "restriction_image_dimension": restriction_rank,
        "supported_code_dimension": len(supported_parameters).bit_length() - 1,
        "supported_code_parameters": supported_parameters,
        "checked_phase_count": phase_count,
        "b_plus_V_all_bent": True,
        "all_partial_spectra_match_signed_translation": True,
        "elapsed_seconds": time.perf_counter() - started,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"Checked {phase_count} phases, period groups, and the four-plane partition.")


if __name__ == "__main__":
    main()
