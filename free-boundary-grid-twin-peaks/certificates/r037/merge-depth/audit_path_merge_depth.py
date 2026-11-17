#!/usr/bin/env python3
"""Exact audit for the path branches used in r051.

No graph enumeration is used.  For P_k the rooted one-peak counts are
binomial(k-1,j), and the sole merge-facing vertex is an endpoint.
The script checks the resulting first and second depth moments by Fraction
summation for 2 <= k <= 40.
"""

from fractions import Fraction
from math import comb


def path_depth_moments(k: int):
    total = sum(comb(k - 1, j) for j in range(k))
    first = sum(Fraction(comb(k - 1, j) * (k - 1 - j), total) for j in range(k))
    second = sum(
        Fraction(comb(k - 1, j) * (k - 1 - j) ** 2, total)
        for j in range(k)
    )
    return total, first, second


def main():
    checked = []
    for k in range(2, 41):
        total, first, second = path_depth_moments(k)
        expected_first = Fraction(k - 1, 2)
        expected_second = Fraction(k * (k - 1), 4)
        assert total == 2 ** (k - 1)
        assert first == expected_first
        assert second == expected_second
        checked.append(k)

    # The two P_k arms in the top-row state have independent roots.
    # E[(k+1 + J_B-J_A)^2] is the stated conditional squared separation.
    for k in checked:
        arm_variance_sum = Fraction(k - 1, 2)
        separation = Fraction((k + 1) ** 2, 1) + arm_variance_sum
        assert separation == Fraction(2 * k * k + 5 * k + 1, 2)

    # A bounded-multiplicity root-to-interface injection cannot exist:
    # endpoint count is 1 whereas some interior root count is at least the
    # average of the k binomial counts.
    for k in checked:
        maximum = max(comb(k - 1, j) for j in range(k))
        assert maximum >= Fraction(2 ** (k - 1), k)

    print("exact Fraction audit passed: P_k endpoint-depth moments, k=2..40")
    print("E[D^2]=k(k-1)/4; two-arm E[rho^2]=(k+1)^2+(k-1)/2")


if __name__ == "__main__":
    main()
