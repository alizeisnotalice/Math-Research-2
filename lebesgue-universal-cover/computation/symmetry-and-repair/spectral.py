"""Certify the homothety of the published degree-five host for spectral restrictions."""
from fractions import Fraction
from pathlib import Path
import json
import sys
from flint import arb, ctx

def exact(value):
    q = Fraction(value)
    return arb(q.numerator) / q.denominator

def run():
    ctx.prec = 224
    area = exact("0.844080666618632325284184")
    budget = exact("0.844082046521218581802948589770767704354")
    def mean_cos(degree):
        return (arb.pi() / (degree + 2)).cos()
    records = []
    for degree in (2455, 2456, 4096):
        c = mean_cos(degree)
        allowed = (c * (budget / area).sqrt() - 1) / 2
        trial_delta = exact("0.00000026") if degree == 4096 else arb(0)
        outer = area * ((1 + 2 * trial_delta) / c) ** 2
        if degree == 2455:
            assert outer > budget
        else:
            assert outer < budget and allowed > trial_delta
        records.append(dict(degree=degree, rotation_mean=str(c),
                            delta_max=str(allowed), trial_delta=str(trial_delta),
                            host_area_at_trial=str(outer), budget_margin=str(budget - outer)))
    return dict(precision_bits=ctx.prec, single_kernel_threshold=2456, degrees=records)

if __name__ == "__main__":
    data = json.dumps(run(), indent=2) + "\n"
    if len(sys.argv) == 2:
        Path(sys.argv[1]).write_text(data)
    else:
        print(data, end="")
