#!/usr/bin/env python3
"""Exact regression for the four weighted signed endpoint budgets.

This checks the identities on the genuine binary product
  A=G=1+x, B=Q=1+x^2, C=AB,
with P=1.  It is a regression of the coefficient bookkeeping only, not a
proof of the all-degree lemma (which is in the accompanying note).
"""
from fractions import Fraction as F


def split(z): return max(z, 0), max(-z, 0)


def truncate_division(num, den, count):
    """First `count` terms of num/den in Q[[x]], den[0]=1."""
    out = []
    for k in range(count):
        nk = F(num[k]) if k < len(num) else F(0)
        out.append(nk - sum(F(den[j]) * out[k-j] for j in range(1, min(k, len(den)-1)+1)))
    return out


def star(poly):
    return list(reversed(poly))


def check(r, b, n, label):
    # Assumes [x^n](R B)=0, r[0]=b[0]=b[n]=1.
    lhs = sum(F(split(r[j])[1]) * b[n-j] for j in range(1, n+1))
    rhs = 1 + sum(F(split(r[j])[0]) * b[n-j] for j in range(1, n+1))
    assert lhs == rhs
    return (label, lhs, rhs)


def conv(p, q):
    out = [F(0)] * (len(p) + len(q) - 1)
    for i, pi in enumerate(p):
        for j, qj in enumerate(q):
            out[i+j] += F(pi) * F(qj)
    return out


def main():
    g, q = [1, 1], [1, 0, 1]
    a, b = g, q
    c = conv(a, b)
    assert c == [1, 1, 1, 1]  # Genuine binary product.
    m, n = len(a)-1, len(b)-1
    r = truncate_division(g, q, n+1)
    s = truncate_division(q, g, m+1)
    rsharp = truncate_division(star(g), star(q), n+1)
    ssharp = truncate_division(star(q), star(g), m+1)
    receipts = [
        check(r, b, n, "R at right endpoint"),
        check(s, a, m, "S at right endpoint"),
        check(rsharp, star(b), n, "R# after simultaneous reversal"),
        check(ssharp, star(a), m, "S# after simultaneous reversal"),
    ]
    for item in receipts:
        print(item)
    # Sharpness: an internal single -1 is possible when the required B entry is 1.
    g2, q2 = [1, 0, 0, 1], [1, 1, 1]
    assert conv(g2, q2) == [1, 1, 1, 1, 1, 1]
    r2 = truncate_division(g2, q2, 3)
    assert r2 == [1, -1, 0]
    sharp = check(r2, q2, 2, "sharpness: internal -1 paid by B_1=1")
    assert sharp[1:] == (F(1), F(1))
    print(sharp)
    print("PASS: all four exact positive/negative endpoint decompositions")


if __name__ == "__main__":
    main()
