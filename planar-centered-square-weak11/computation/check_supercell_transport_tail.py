#!/usr/bin/env python3
"""Exact safety check for the supercell transport tail bound.

It reads a frozen periodic input and checks the bound proved in
proofs/supercell-transport-tail.md.  The program deliberately does not build
or optimise a supercell; it only certifies the range in which any labelled
atom matching with l-infinity displacement at most delta retains the r < R
complete-halo reduction.
"""
from __future__ import annotations

import argparse
from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys


def fraction(value: str) -> Fraction:
    return Fraction(value)


def bound(mass: Fraction, px: Fraction, py: Fraction, radius: Fraction,
          delta: Fraction) -> Fraction:
    """rho(1 + (Px+2d)/(2r))(1 + (Py+2d)/(2r))."""
    rho = mass / (px * py)
    return rho * (1 + (px + 2 * delta) / (2 * radius)) * (
        1 + (py + 2 * delta) / (2 * radius)
    )


def ceil_fraction(value: Fraction) -> int:
    return -((-value.numerator) // value.denominator)


def check_nonunit_lambda_regression() -> dict[str, str]:
    """Keep the lambda factor in the diagnostic safe-root formula."""
    mass = px = py = radius = Fraction(1)
    lam = Fraction(3)
    delta = Fraction(1, 5)
    rho = mass / (px * py)
    a, b = 2 * radius + px, 2 * radius + py
    c = 4 * radius * radius * lam / rho
    assert bound(mass, px, py, radius, delta) < lam
    assert (a + 2 * delta) * (b + 2 * delta) < c
    assert c == 12  # Omitting lam would incorrectly give 4.
    return {
        "lambda": str(lam),
        "delta": str(delta),
        "C": str(c),
        "tail": str(bound(mass, px, py, radius, delta)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        default="planar-centered-square-weak11/computation/inputs/r004_final_0504.json",
    )
    parser.add_argument("--radius", default="5")
    parser.add_argument("--safe-delta", default="1/32")
    parser.add_argument("--reject-delta", default="1/31")
    args = parser.parse_args()
    regression = check_nonunit_lambda_regression()

    path = Path(args.input)
    raw = path.read_bytes()
    data = json.loads(raw)
    px, py, lam = (fraction(data[key]) for key in ("Px", "Py", "lam"))
    mass = sum(fraction(point[2]) for point in data["pts"])
    radius = fraction(args.radius)
    safe_delta = fraction(args.safe_delta)
    reject_delta = fraction(args.reject_delta)
    assert radius > 0 and safe_delta >= 0 and reject_delta >= safe_delta
    safe = bound(mass, px, py, radius, safe_delta)
    rejected = bound(mass, px, py, radius, reject_delta)
    baseline = bound(mass, px, py, radius, Fraction(0))
    assert baseline < lam
    assert safe < lam
    assert rejected >= lam
    rho = mass / (px * py)
    # Given baseline < lambda, (A+2d)(B+2d)<C describes the exact
    # positive safe interval for any positive lambda.
    a, b, c = 2 * radius + px, 2 * radius + py, 4 * radius * radius * lam / rho
    discriminant = (a - b) ** 2 + 4 * c
    with localcontext() as context:
        context.prec = 70
        decimal_delta_star = (
            - (Decimal(a.numerator) / Decimal(a.denominator)
               + Decimal(b.numerator) / Decimal(b.denominator))
            + (Decimal(discriminant.numerator)
               / Decimal(discriminant.denominator)).sqrt()
        ) / 4
    # The formula decreases in r because both positive factors do, so this
    # one exact endpoint comparison certifies every r >= radius.
    result = {
        "status": "passed",
        "input": str(path),
        "input_sha256": hashlib.sha256(raw).hexdigest(),
        "python": sys.version.split()[0],
        "atoms": len(data["pts"]),
        "lambda": str(lam),
        "period_x": str(px),
        "period_y": str(py),
        "cell_mass": str(mass),
        "radius": str(radius),
        "safe_delta": str(safe_delta),
        "baseline_tail_at_zero_displacement": str(baseline),
        "safe_interval": {
            "condition": "(A+2*delta)*(B+2*delta) < C",
            "A": str(a),
            "B": str(b),
            "C": str(c),
            "delta_star_decimal_60": format(decimal_delta_star, ".60f"),
        },
        "transport_tail_at_safe_delta": str(safe),
        "safe_slack": str(lam - safe),
        "reject_delta": str(reject_delta),
        "transport_tail_at_reject_delta": str(rejected),
        "reject_excess": str(rejected - lam),
        "supercell_q": ceil_fraction(radius / min(2 * px, py)) + 1,
        "scope": (
            "For every labelled 2Px-by-Py supercell perturbation whose "
            "atom matching moves each mass by l-infinity distance at most "
            "safe_delta, no strict lambda witness has radius >= radius. "
            "The reject_delta result only shows that this particular "
            "uniform transport bound no longer proves the claim."
        ),
        "nonunit_lambda_regression": regression,
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
