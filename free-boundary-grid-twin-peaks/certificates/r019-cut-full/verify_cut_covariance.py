#!/usr/bin/env python3
"""Exact audit of r019 cut-covariance JSON files."""
import json
from fractions import Fraction
from pathlib import Path


BASE = Path(__file__).parent
ROOT_MOMENTS = BASE.parent / "r018-root-moments"
INDEPENDENT_CUT = BASE.parent / "r019-association"
EXPECTED_A2 = {
    4: 779299174080,
    5: 7063261481362958280000,
    6: 304351937651819767235871772704107520,
}


def load(name):
    return json.loads((BASE / name).read_text())


def main():
    ie3 = load("n3-ie.json")
    direct3 = load("n3-permutation-audit.json")
    for key in ("A2", "cut_marginal_numerators", "ordered_pair_cross_numerators",
                "covariance_common_denominator", "covariance_numerators"):
        assert ie3[key] == direct3[key], key
    summary = [{"n": 3, "verification": "all 9! label permutations equals IE"}]
    for n, a2_expected in EXPECTED_A2.items():
        data = load(f"n{n}-cut-covariance.json")
        a2 = int(data["A2"])
        assert a2 == a2_expected
        L = list(map(int, data["cut_marginal_numerators"]))
        H = [[int(v) for v in row] for row in data["ordered_pair_cross_numerators"]]
        C = [[int(v) for v in row] for row in data["covariance_numerators"]]
        denominator = int(data["covariance_common_denominator"])
        assert denominator == 4 * a2 * a2
        assert len(L) == len(H) == len(C) == n - 1
        assert H == [list(row) for row in zip(*H)]
        assert C == [list(row) for row in zip(*C)]
        for j in range(n - 1):
            for k in range(n - 1):
                assert C[j][k] == 2 * a2 * H[j][k] - L[j] * L[k]
        if n <= 5:
            # Independently written r019 accumulator, with different field names.
            other = json.loads((INDEPENDENT_CUT / f"n{n}-cut-ie.json").read_text())
            assert int(other["A2"]) == a2
            assert list(map(int, other["cut_mean_numerators"])) == L
            assert [[int(v) for v in row]
                    for row in other["cut_ordered_pair_numerators"]] == H

        # sum_j s_j(x)=2x-(n-1).  The matrix sum is therefore the covariance
        # of doubled x coordinates. D4 gives equal x/y covariances, hence
        # C_dot=G/(4A2)=matrix_sum/2 and one original-coordinate term is
        # matrix_sum/4.
        root = json.loads((ROOT_MOMENTS / f"n{n}-root-moments.json").read_text())
        assert int(root["A2"]) == a2
        G = int(root["secondsum"])
        matrix_sum = Fraction(sum(sum(row) for row in C), denominator)
        c_dot = Fraction(G, 4 * a2)
        c_x = matrix_sum / 4
        assert c_dot == matrix_sum / 2
        assert c_x == c_dot / 2
        positives = [{"j": j, "k": k, "numerator": str(C[j][k]),
                      "decimal": float(Fraction(C[j][k], denominator))}
                     for j in range(n - 1) for k in range(n - 1) if C[j][k] > 0]
        summary.append({"n": n, "A2": str(a2), "all_entries_nonpositive": not positives,
                        "positive_entries": positives,
                        "means": [str(Fraction(x, 2 * a2)) for x in L],
                        "independent_n3_to_n5_cut_accumulator_agrees": True if n <= 5 else None,
                        "matrix_sum": str(matrix_sum),
                        "original_coordinate_covariance": str(c_x),
                        "r018_C_dot": str(c_dot)})
    output = {"status": "passed", "runs": summary}
    (BASE / "cut-covariance-verified.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
