#!/usr/bin/env python3
"""Exact 3x3 enumeration: a minimal counterexample to peak-polynomial real-rootedness.

The vertices are labelled by the entries of a permutation of 0,...,8.  A peak
is a vertex larger than every grid-neighbour.  This deliberately does not use
the independent-set recurrence, so it independently verifies the displayed
histogram used in the polynomial obstruction.
"""
from itertools import permutations
from math import comb, factorial


def neighbours(n):
    ans = []
    for r in range(n):
        for c in range(n):
            q = []
            for dr, dc in ((-1, 0), (1, 0), (0, -1), (0, 1)):
                rr, cc = r + dr, c + dc
                if 0 <= rr < n and 0 <= cc < n:
                    q.append(rr * n + cc)
            ans.append(tuple(q))
    return ans


def main():
    n = 3
    adj = neighbours(n)
    hist = [0] * (n * n + 1)
    for lab in permutations(range(n * n)):
        k = sum(all(lab[v] > lab[w] for w in adj[v]) for v in range(n * n))
        hist[k] += 1
    expected = [0, 29568, 143808, 160704, 24000, 4800, 0, 0, 0, 0]
    assert hist == expected, (hist, expected)
    assert sum(hist) == factorial(9)

    # Q(z)=F(z)/z has degree four and coefficients q_0,...,q_4.
    q = hist[1:6]
    assert q[3] * q[3] < q[2] * q[4]  # ordinary log-concavity already fails
    b2 = q[2] // comb(4, 2)
    b3 = q[3] // comb(4, 3)
    b4 = q[4] // comb(4, 4)
    assert b3 * b3 < b2 * b4  # Newton's necessary inequality
    print("histogram", hist[1:6])
    print("ordinary_LC_failure", q[3] * q[3], "<", q[2] * q[4])
    print("Newton_failure", b3 * b3, "<", b2 * b4)
    print("conclusion: F_3(z)/z is not real-rooted; neither is F_3(z).")


if __name__ == "__main__":
    main()
