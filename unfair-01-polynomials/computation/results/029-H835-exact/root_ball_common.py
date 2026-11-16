"""Exact rational helpers for the R8--L9 root-ball certificate prototype.

No floating-point arithmetic is used in this module.  Complex numbers are
pairs of fractions and every absolute-value bound is produced from integer
square roots at a declared rational scale.
"""
from __future__ import annotations

from fractions import Fraction
from math import comb, isqrt


ComplexQ = tuple[Fraction, Fraction]


def qcomplex(re: Fraction | int = 0, im: Fraction | int = 0) -> ComplexQ:
    return Fraction(re), Fraction(im)


def cadd(a: ComplexQ, b: ComplexQ) -> ComplexQ:
    return a[0] + b[0], a[1] + b[1]


def csub(a: ComplexQ, b: ComplexQ) -> ComplexQ:
    return a[0] - b[0], a[1] - b[1]


def cmul(a: ComplexQ, b: ComplexQ) -> ComplexQ:
    return a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]


def cscale(a: ComplexQ, q: Fraction | int) -> ComplexQ:
    return a[0] * q, a[1] * q


def cpow(a: ComplexQ, n: int) -> ComplexQ:
    out = qcomplex(1, 0)
    base = a
    while n:
        if n & 1:
            out = cmul(out, base)
        base = cmul(base, base)
        n //= 2
    return out


def csqabs(a: ComplexQ) -> Fraction:
    return a[0] * a[0] + a[1] * a[1]


def fraction_string(q: Fraction) -> str:
    return f"{q.numerator}/{q.denominator}"


def parse_fraction(text: str) -> Fraction:
    numerator, denominator = text.split("/", 1)
    return Fraction(int(numerator), int(denominator))


def complex_from_strings(pair: list[str]) -> ComplexQ:
    return Fraction(pair[0]), Fraction(pair[1])


def complex_to_json(a: ComplexQ) -> dict[str, str]:
    return {"re": fraction_string(a[0]), "im": fraction_string(a[1])}


def complex_from_json(value: dict[str, str]) -> ComplexQ:
    return parse_fraction(value["re"]), parse_fraction(value["im"])


def abs_bounds(a: ComplexQ, scale: int) -> tuple[Fraction, Fraction]:
    """Return exact rational lower/upper bounds for |a|.

    If t=floor(scale^2 |a|^2), then isqrt(t)/scale is a lower
    bound and (isqrt(t)+1)/scale is a strict upper bound.  This uses
    only integer operations on the rational squared modulus.
    """
    if scale <= 0:
        raise ValueError("positive scale required")
    squared = csqabs(a)
    if squared < 0:
        raise AssertionError("squared modulus cannot be negative")
    floor_scaled_squared = (squared.numerator * scale * scale) // squared.denominator
    lower_integer = isqrt(floor_scaled_squared)
    return Fraction(lower_integer, scale), Fraction(lower_integer + 1, scale)


def taylor_coefficients(coefficients: list[int], center: ComplexQ) -> list[ComplexQ]:
    """Taylor coefficients f^(k)(center)/k! for an ascending polynomial."""
    answer: list[ComplexQ] = []
    for k in range(len(coefficients)):
        total = qcomplex()
        for n in range(k, len(coefficients)):
            if coefficients[n]:
                total = cadd(total, cscale(cpow(center, n - k), coefficients[n] * comb(n, k)))
        answer.append(total)
    return answer


def polynomial_rouche_bounds(coefficients: list[int], center: ComplexQ, radius: Fraction,
                              scale: int) -> tuple[Fraction, Fraction]:
    """Bounds for the linear Rouché test at one disk.

    Returns LHS upper and RHS lower for
      |f(c)| + sum_{k>=2}|f^(k)(c)/k!|r^k < |f'(c)|r.
    """
    terms = taylor_coefficients(coefficients, center)
    _, lhs = abs_bounds(terms[0], scale)
    for k in range(2, len(terms)):
        _, upper = abs_bounds(terms[k], scale)
        lhs += upper * radius**k
    lower, _ = abs_bounds(terms[1], scale)
    rhs = lower * radius
    return lhs, rhs


def lift_rouche_bounds(alpha_center: ComplexQ, alpha_radius: Fraction, center: ComplexQ,
                       radius: Fraction, scale: int, exponent: int) -> tuple[Fraction, Fraction]:
    """Rouché bounds for z^exponent-alpha with alpha in a mother disk."""
    if exponent < 2:
        raise ValueError("lift exponent must be at least two")
    g0 = csub(cpow(center, exponent), alpha_center)
    _, lhs = abs_bounds(g0, scale)
    lhs += alpha_radius
    # z^q Taylor coefficient at order k is binom(q,k)c^(q-k).
    for k in range(2, exponent + 1):
        _, upper = abs_bounds(cscale(cpow(center, exponent - k), comb(exponent, k)), scale)
        lhs += upper * radius**k
    lower, _ = abs_bounds(cscale(cpow(center, exponent - 1), exponent), scale)
    rhs = lower * radius
    return lhs, rhs


def disk_disjoint(a: ComplexQ, ra: Fraction, b: ComplexQ, rb: Fraction) -> bool:
    """Exact squared-distance test for open disks."""
    return csqabs(csub(a, b)) > (ra + rb) ** 2


def disk_nonreal(center: ComplexQ, radius: Fraction) -> bool:
    return abs(center[1]) > radius
