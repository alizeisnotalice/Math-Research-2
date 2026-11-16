#!/usr/bin/env python3
"""Independent exact checks for D4-compressed single-root counts."""
from fractions import Fraction
import hashlib
import json
from pathlib import Path


BASE = Path(__file__).parent
SMALL = BASE / "reference"
KNOWN_N7 = json.loads((SMALL / "n7-d4-root-moments.json").read_text())


def at(roots, n, x, y):
    return roots[y * n + x]


def audit(n):
    data = json.loads((BASE / f"n{n}.json").read_text())
    roots = list(map(int, data["root_counts"]))
    assert len(roots) == n * n
    assert sum(roots) == int(data["A1"])
    if n in (3, 4):
        reference = json.loads((SMALL / f"n{n}-roots.json").read_text())
        assert roots == list(map(int, reference["single_peak_root_counts"]))
    if n == 7:
        for key in ("A1", "A2", "squared_euclidean_total"):
            assert data[key] == KNOWN_N7[key]
    for y in range(n):
        for x in range(n):
            assert at(roots, n, x, y) == at(roots, n, n - 1 - x, y)
            assert at(roots, n, x, y) == at(roots, n, x, n - 1 - y)
    coordinate_edges = 0
    for y in range(n):
        for x in range((n - 1) // 2):
            assert at(roots, n, x + 1, y) >= at(roots, n, x, y)
            coordinate_edges += 1
    for x in range(n):
        for y in range((n - 1) // 2):
            assert at(roots, n, x, y + 1) >= at(roots, n, x, y)
            coordinate_edges += 1
    all_row = max(
        Fraction(at(roots, n, x + 1, y), at(roots, n, x, y))
        * Fraction((x + 1) ** 2, n * n)
        for y in range(n) for x in range(n - 1)
    )
    assert all_row <= 1
    strict_binomial = max(
        (
            Fraction(at(roots, n, x + 1, y), at(roots, n, x, y))
            / Fraction((n - 1 - x) ** 2, (x + 1) ** 2)
            for y in range(n) for x in range((n - 1) // 2)
        ),
        default=Fraction(0),
    )
    assert strict_binomial <= 1
    central = 0
    if n % 2 == 0:
        x = n // 2 - 1
        for y in range(n):
            assert at(roots, n, x, y) == at(roots, n, x + 1, y)
            central += 1
    return {
        "n": n,
        "A1": data["A1"],
        "A2": data["A2"],
        "root_sum_matches_A1": True,
        "small_full_vector_match": n in (3, 4),
        "n7_historical_A1_A2_D_match": n == 7,
        "strict_coordinate_edges_checked": coordinate_edges,
        "all_row_N_over_i_plus_1_squared_constant": str(all_row),
        "strict_binomial_square_constant": str(strict_binomial),
        "central_reflection_edges_checked": central,
    }


def main():
    result = {
        "status": "passed",
        "source_sha256": hashlib.sha256(
            (BASE / "root_counts_d4.cpp").read_bytes()).hexdigest(),
        "runs": [audit(n) for n in (3, 4, 7)],
    }
    (BASE / "root-counts-verified.json").write_text(
        json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
