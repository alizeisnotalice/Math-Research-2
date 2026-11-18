"""Outward-rounded degree-grid losses and least full-turn sample counts."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import platform

import flint
from flint import arb, ctx


def rational(q):
    q = Fraction(q)
    return arb(q.numerator) / q.denominator


def true_loss(distance, tangential_radius):
    # cos(distance/2)-1 in its stable half-angle form.
    return 2 * tangential_radius * (distance / 2).sin() - 2 * (distance / 4).sin() ** 2


def old_point_loss(distance, tangential_radius):
    return (distance / 2).sin() ** 2 + tangential_radius * distance.sin()


def decimal_upper(value, digits):
    scale = 10**digits
    z = (value.upper() * scale).ceil().unique_fmpz()
    assert z is not None
    return f"{int(z)//scale}.{int(z)%scale:0{digits}d}"


def main(output):
    assert __debug__
    ctx.prec = 256
    pi = arb.pi()
    root3 = arb(3).sqrt()
    radius = 1 / root3
    slope = 1 / (2 * root3)
    grids = []
    for k in range(13):
        spacing = Fraction(1, 10**k)
        distance = rational(spacing) * pi / 360
        truth = true_loss(distance, slope)
        material = old_point_loss(distance, slope)
        new_plane = 2 * (distance / 2).sin() ** 2
        assert 0 < new_plane < truth < material
        grids.append(
            dict(
                mesh_degrees=str(spacing),
                nearest_radians=str(distance),
                true_support_loss=str(truth),
                old_material_loss=str(material),
                new_support_minus_old_point=str(new_plane),
                upper_decimals=[
                    decimal_upper(x, 32) for x in (truth, material, new_plane)
                ],
            )
        )
    resolutions = []
    for k in range(1, 16):
        eps = Fraction(1, 10**k)
        tol = rational(eps)
        critical_distance = 2 * (pi / 6 - ((1 + tol) / (2 * radius)).acos())
        n = (pi / critical_distance).floor().unique_fmpz()
        assert n is not None
        count = int(n) + 1
        assert pi / (count - 1) < pi / 3
        assert true_loss(pi / count, slope) < tol < true_loss(pi / (count - 1), slope)
        resolutions.append(
            dict(
                tolerance=str(eps),
                max_mesh_degrees=str((360 / pi) * critical_distance),
                strict_min_full_turn_count=count,
            )
        )
    results = dict(
        precision_bits=ctx.prec,
        python=platform.python_version(),
        python_flint=flint.__version__,
        sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        grids=grids,
        resolutions=resolutions,
    )
    output.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    main(parser.parse_args().output)
