"""Directed bounds for the regular-normal support obstruction and its angular tolerance."""
from fractions import Fraction
from pathlib import Path
import json
import sys

from flint import arb, ctx

ctx.prec = 192


def fraction_ball(value):
    value = Fraction(value)
    return arb(value.numerator) / value.denominator


def floor_decimal(value, places=20):
    factor = 10**places
    scaled = (value.lower() * factor).floor().unique_fmpz()
    result = int(scaled)
    return f"{result // factor}.{result % factor:0{places}d}"


amplitude = fraction_ball("1/16")
records = []
for degrees in ("0", ".001", ".1", "1", "2"):
    alpha = fraction_ball(degrees) * arb.pi() / 180
    angular_cover = (arb.pi()/5 + alpha).cos()
    assert angular_cover > 0
    scaled_error = alpha / angular_cover
    fraction = (fraction_ball("1/2") - 3*alpha - scaled_error) / (1 + scaled_error)
    release = amplitude * fraction
    assert release > 0
    records.append({"angular_tolerance_degrees": degrees,
                    "release_lower_floor": floor_decimal(release),
                    "required_maximum_support_floor": floor_decimal(fraction_ball("1/2") + release)})
# The sharp pentagon pose has four active supports and a slack fifth.
for i in range(1, 5):
    angle = 2 * arb.pi() * i / 5
    assert (angle.cos() + (3*angle).cos() + fraction_ball("1/2")).contains(0)
assert -2*amplitude < amplitude/2
result = {"precision_bits": ctx.prec, "cubic_amplitude": "1/16",
          "regular_pentagon_exact_uniform_release": "1/32",
          "regular_pentagon_support_minimum": "17/32", "perturbations": records}
Path(sys.argv[1]).write_text(json.dumps(result, indent=2) + "\n")
print(records)
