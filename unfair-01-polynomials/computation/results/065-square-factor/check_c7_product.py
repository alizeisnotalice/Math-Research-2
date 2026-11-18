"""Apply proved local necessity tests to a binary product in the C7-gap branch."""

import argparse
import json
from math import comb
from pathlib import Path

import sympy as sp


Z = sp.Symbol("z")
CONSTANTS = {
    163: (124, [1, 6, 121, 106, 68, 115, 117]),
    1279: (1133, [1, 1111, 1202, 215, 1193, 762, 576]),
}


def hasse(poly, order):
    return sp.Poly.from_dict(
        {(exponent - order,): coeff * comb(exponent, order)
         for (exponent,), coeff in poly.terms() if exponent >= order},
        (Z,), domain=poly.domain,
    )


def over_integers(poly):
    return sp.Poly.from_dict(poly.as_dict(), (Z,), domain=sp.ZZ)


def power_mod(value, exponent, modulus):
    result = sp.Poly(1, Z, domain=modulus.domain)
    while exponent:
        if exponent & 1:
            result = (result * value).rem(modulus)
        value = (value * value).rem(modulus)
        exponent //= 2
    return result


def verify_local(g, prime, squarefree_product):
    h, prefix = CONSTANTS[prime]
    residue = g.set_modulus(prime)
    _, factors = residue.factor_list()
    square_degree = sum(factor.degree() * (count // 2) for factor, count in factors)
    report = dict(prime=prime, square_degree=square_degree)
    if square_degree < 7:
        return dict(report, status="excluded", reason="square_degree")
    if square_degree > 7:
        return dict(report, status="unresolved", reason="above_boundary")
    j = sp.Poly(1, Z, modulus=prime)
    for factor, count in factors:
        j *= factor ** (count // 2)
    j_coeff = [int(coeff) % prime for coeff in j.all_coeffs()]
    report["j_descending"] = j_coeff
    if j_coeff[:7] != prefix or not j_coeff[-1]:
        return dict(report, status="excluded", reason="j_prefix")

    checks = []
    for factor, count in factors:
        mu = count // 2
        if not mu:
            continue
        denominator = hasse(j, mu).rem(factor)
        assert not denominator.is_zero
        beta = (sp.Poly(Z**6, Z, modulus=prime) * denominator.invert(factor)).rem(factor)
        lhs = (h * beta**2 * hasse(residue, count)).rem(factor)
        assert not lhs.is_zero
        # Reduction by an integral lift of the monic factor evaluates at
        # an unramified lift of any one of its roots.
        value = hasse(g, count - 2 * mu).rem(over_integers(factor))
        assert all(int(coeff) % prime == 0 for coeff in value.coeffs())
        rhs = sp.Poly.from_dict(
            {key: int(coeff) // prime for key, coeff in value.as_dict().items()},
            (Z,), modulus=prime,
        )
        if squarefree_product:
            passes = lhs == rhs
        elif rhs.is_zero:
            passes = False
        else:
            ratio = (lhs * rhs.invert(factor)).rem(factor)
            passes = ratio.degree() == 0 and power_mod(ratio, (prime - 1) // 2, factor).is_one
        checks.append(dict(
            factor_descending=[int(coeff) % prime for coeff in factor.all_coeffs()],
            multiplicity=count, passes=passes,
        ))
    return dict(report, status="unresolved" if all(row["passes"] for row in checks) else "excluded",
                reason="lift", checks=checks)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="JSON list of binary coefficients in ascending order")
    args = parser.parse_args()
    values = json.loads(args.input.read_text())
    if (not isinstance(values, list) or len(values) < 12
            or any(type(value) is not int or value not in (0, 1) for value in values)
            or values[:12] != [1] * 7 + [0] + [1] * 4 or values[-1] != 1):
        parser.error("input requires the normalized C7-gap prefix and leading coefficient 1")
    c = sp.Poly.from_list(values[::-1], Z, domain=sp.ZZ)
    f = c.sqf_part()
    constant = int(f.nth(0))
    assert constant in (-1, 1)
    g = sp.Poly.from_list([constant * (-1)**i * f.nth(i) for i in range(f.degree() + 1)], Z, domain=sp.ZZ)
    assert g.is_monic
    print(json.dumps([verify_local(g, p, f == c) for p in CONSTANTS], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
