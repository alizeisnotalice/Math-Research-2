#!/usr/bin/env python3
"""Exact tail check for the mass-dominated supercell transport lemma.

This checker verifies only the uniform tail prerequisite.  A subsequent halo
calculation must still use the target supercell's actual total mass and its
complete rectangle union.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys


def q(value: str) -> Fraction:
    return Fraction(value)


def unscaled_tail(mass: Fraction, px: Fraction, py: Fraction,
                  radius: Fraction, delta: Fraction) -> Fraction:
    rho = mass / (px * py)
    return rho * (1 + (px + 2 * delta) / (2 * radius)) * (
        1 + (py + 2 * delta) / (2 * radius)
    )


def assert_positive_parameters(px: Fraction, py: Fraction, lam: Fraction) -> None:
    assert px > 0
    assert py > 0
    assert lam > 0


def check_bad_parameter_regression() -> dict[str, list[str]]:
    """Exercise the same positivity assertions used by the command entrypoint."""
    rejected = []
    for label, px, py, lam in (
        ("zero_period_x", Fraction(0), Fraction(1), Fraction(1)),
        ("negative_period_y", Fraction(1), Fraction(-1), Fraction(1)),
        ("zero_lambda", Fraction(1), Fraction(1), Fraction(0)),
    ):
        try:
            assert_positive_parameters(px, py, lam)
        except AssertionError:
            rejected.append(label)
        else:
            raise AssertionError(f"invalid parameters accepted: {label}")
    assert len(rejected) == 3
    return {"rejected": rejected}


def check_nonunit_lambda_regression() -> dict[str, str]:
    """A purely rational regression for the factor kappa and lambda."""
    mass = px = py = radius = Fraction(1)
    delta, kappa, lam = Fraction(1, 5), Fraction(7, 8), Fraction(3)
    base = unscaled_tail(mass, px, py, radius, delta)
    assert base == Fraction(289, 100)
    assert kappa * base < lam
    assert kappa < lam / base
    return {
        "lambda": str(lam), "kappa": str(kappa), "delta": str(delta),
        "unscaled_tail": str(base), "scaled_tail": str(kappa * base),
        "kappa_max": str(lam / base),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        default="planar-centered-square-weak11/computation/inputs/r008_final_0611.json",
    )
    parser.add_argument("--radius", default="5")
    parser.add_argument("--kappa", default="1001/1000")
    parser.add_argument("--delta", action="append", default=[])
    parser.add_argument("--self-test-bad-parameters", action="store_true")
    args = parser.parse_args()
    bad_parameter_regression = check_bad_parameter_regression()
    if args.self_test_bad_parameters:
        print(json.dumps({"status": "passed", "bad_parameter_regression": bad_parameter_regression},
                         indent=2, sort_keys=True))
        return
    deltas = [q(value) for value in (args.delta or ["1/1000", "1/10000"])]

    path = Path(args.input)
    raw = path.read_bytes()
    data = json.loads(raw)
    px, py, lam = (q(data[key]) for key in ("Px", "Py", "lam"))
    assert_positive_parameters(px, py, lam)
    weights = [q(point[2]) for point in data["pts"]]
    mass, radius, kappa = sum(weights), q(args.radius), q(args.kappa)
    assert mass > 0 and radius > 0 and kappa > 0 and all(weight > 0 for weight in weights)
    checks = []
    for delta in deltas:
        assert delta >= 0
        base = unscaled_tail(mass, px, py, radius, delta)
        scaled = kappa * base
        kappa_max = lam / base
        assert scaled < lam
        assert kappa < kappa_max
        checks.append({
            "delta": str(delta),
            "G_delta_at_radius": str(base),
            "kappa_max_exclusive": str(kappa_max),
            "kappa_times_G": str(scaled),
            "strict_slack": str(lam - scaled),
        })
    print(json.dumps({
        "status": "passed",
        "input": str(path),
        "input_sha256": hashlib.sha256(raw).hexdigest(),
        "python": sys.version.split()[0],
        "atoms": len(weights),
        "period_x": str(px), "period_y": str(py), "lambda": str(lam),
        "base_cell_mass": str(mass), "radius": str(radius), "kappa": str(kappa),
        "checks": checks,
        "bad_parameter_regression": bad_parameter_regression,
        "nonunit_lambda_regression": check_nonunit_lambda_regression(),
        "scope": (
            "Every labelled target atom with mass at most kappa times its "
            "positive base-label mass and l-infinity displacement at most "
            "delta has no strict witness of radius >= radius. This does not "
            "compute the target supercell's halo or its actual total mass."
        ),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
