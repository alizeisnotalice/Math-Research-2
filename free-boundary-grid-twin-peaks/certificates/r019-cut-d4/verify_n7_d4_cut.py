#!/usr/bin/env python3
"""Exact audit for D4-orbit cut covariance, including n=7."""
from fractions import Fraction
import json
from pathlib import Path


BASE = Path(__file__).parent
FULL = BASE.parent / "r019-cut-full"
ROOT = BASE.parent / "r018-root-moments"
EXPECTED_A2 = {
    4: 779299174080,
    5: 7063261481362958280000,
    6: 304351937651819767235871772704107520,
    7: 137077928808481869759712347602605837601908319657172480,
}


def load(path):
    return json.loads(path.read_text())


def main():
    runs = []
    for n, expected_a2 in EXPECTED_A2.items():
        data = load(BASE / f"n{n}-d4-cut.json")
        a2 = int(data["A2"])
        assert a2 == expected_a2
        L = list(map(int, data["cut_marginal_numerators"]))
        H = [[int(x) for x in row] for row in data["ordered_pair_cross_numerators"]]
        C = [[int(x) for x in row] for row in data["covariance_numerators"]]
        den = int(data["covariance_common_denominator"])
        assert den == 4 * a2 * a2
        assert len(L) == len(H) == len(C) == n - 1
        assert H == [list(row) for row in zip(*H)]
        assert C == [list(row) for row in zip(*C)]
        for j in range(n - 1):
            for k in range(n - 1):
                assert C[j][k] == 2 * a2 * H[j][k] - L[j] * L[k]
                assert C[j][k] == C[n - 2 - j][n - 2 - k]

        if n <= 6:
            full = load(FULL / f"n{n}-cut-covariance.json")
            for key in ("A2", "cut_marginal_numerators",
                        "ordered_pair_cross_numerators",
                        "covariance_common_denominator",
                        "covariance_numerators",
                        "all_covariance_entries_nonpositive"):
                assert data[key] == full[key], (n, key)
        else:
            assert data["raw_independent_sets"] == 1280128950

        root_name = f"n{n}-{'d4-' if n == 7 else ''}root-moments.json"
        root = load(ROOT / root_name)
        assert int(root["A2"]) == a2
        matrix_sum = Fraction(sum(sum(row) for row in C), den)
        c_dot = Fraction(int(root["secondsum"]), 4 * a2)
        assert c_dot == matrix_sum / 2
        positives = [{"j": j, "k": k,
                      "numerator": str(C[j][k]),
                      "decimal": float(Fraction(C[j][k], den))}
                     for j in range(n - 1) for k in range(n - 1)
                     if C[j][k] > 0]
        if n == 7:
            assert {(x["j"], x["k"]) for x in positives} == {(0, 5), (5, 0)}
            assert all(C[j][k] < 0 for j in range(n - 1)
                       for k in range(n - 1) if (j, k) not in {(0, 5), (5, 0)})
        runs.append({
            "n": n,
            "A2": str(a2),
            "full_table_matrix_agreement": True if n <= 6 else None,
            "matrix_sum": str(matrix_sum),
            "r018_C_dot": str(c_dot),
            "positive_entries": positives,
        })
    out = {"status": "passed", "runs": runs}
    (BASE / "d4-cut-covariance-verified.json").write_text(
        json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
