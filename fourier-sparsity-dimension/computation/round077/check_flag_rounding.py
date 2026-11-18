"""Reproduce with python3 -B check_flag_rounding.py. Uses only the standard library."""
import json
import math
import platform
from fractions import Fraction as Q


def check_cyclic_columns(n, h):
    mask = (1 << h) - 1
    column = 1
    seen = []
    for _ in range(n):
        seen.append(column)
        column = (column ^ (column << 1)) & mask
    assert column == 1 and len(set(seen)) == n
    assert all(seen[i].bit_length() == i + 1 for i in range(h))


def exact_squared_ratio(k, n, copies_scale):
    # Float only for choosing nearby integer ranks. The reported rational
    # maximum for the resulting explicitly listed instance is exact.
    ranks = [round(n * 2 ** (-1 + (j + 0.5) / k)) for j in range(k)]
    targets = [Q(h, n) for h in ranks]
    copies = [round(copies_scale / (t*t)) for t in targets]
    assert all(Q(1, 2) < t < 1 for t in targets)
    beta_squared = sum(w * t*t for t, w in zip(targets, copies))
    # A common multiplicity n (the number of columns per block) cancels.
    powers = [Q(1, 4), Q(1, 2), Q(1)]
    events = {Q(1), Q(2)}
    for t in targets:
        for low, high in zip(powers, powers[1:]):
            event = 2*t/(low+high)
            if 1 < event < 2:
                events.add(event)
    events = sorted(events)
    best = Q(0)
    for left, right in zip(events, events[1:]):
        mid = (left+right)/2
        selected = [min(powers, key=lambda z: abs(t-mid*z)) for t in targets]
        linear = sum(w*t*z for w, t, z in zip(copies, targets, selected))
        quadratic = sum(w*z*z for w, z in zip(copies, selected))
        scale = min(right, max(left, linear/quadratic))
        best = max(best, 2*linear*scale-quadratic*scale*scale)
    ratio = best/beta_squared
    assert ratio <= 1
    bound = 1/(3*k*math.tanh(math.log(2)/(2*k)))
    assert float(ratio) >= bound - 2e-15
    return {
        "k": k, "n": n, "copies_scale": copies_scale, "ranks": ranks,
        "squared_ratio_exact": str(ratio), "squared_ratio_decimal": float(ratio),
        "finite_phase_bound_decimal": bound
    }


def walsh(values):
    width = 1
    while width < len(values):
        for start in range(0, len(values), 2*width):
            for j in range(start, start+width):
                a, b = values[j], values[j+width]
                values[j], values[j+width] = a+b, a-b
        width *= 2
    return values


def check_boolean_block_spectrum(ranks):
    leaves = len(ranks)
    p = leaves.bit_length()-1
    assert leaves == 2**p
    shifts = [p + sum(ranks[:i]) for i in range(leaves)]
    dimension = p + sum(ranks)
    spectrum = set()
    for h, shift in zip(ranks, shifts):
        for private in range(1, 2**h, 2):
            for address in range(leaves):
                spectrum.add(address + (private << shift))
    truth = []
    for physical in range(2**dimension):
        leaf = physical % leaves
        private = (physical >> shifts[leaf]) % (2**ranks[leaf])
        value = (-1 if private % 2 else 1) * (-1 if private < 2 else 1)
        truth.append(value)
    assert {i for i, coefficient in enumerate(walsh(truth)) if coefficient} == spectrum
    transforms = walsh([int(i in spectrum) for i in range(2**dimension)])
    minimum_tail = (len(spectrum)-max(transforms[1:]))//2
    expected = leaves * 2**(min(ranks)-2)
    assert minimum_tail == expected
    return {"private_ranks": ranks, "effective_rank": dimension,
            "sparsity": len(spectrum), "minimum_hyperplane_tail": minimum_tail}


def main():
    cases = [(n, h) for n in (8, 16, 32) for h in range(n//2+1, n+1)]
    for n, h in cases:
        check_cyclic_columns(n, h)
    print(json.dumps({
        "python": platform.python_version(), "cyclic_column_cases": len(cases),
        "boolean_cases": [check_boolean_block_spectrum(h) for h in ([3, 4], [3, 3, 4, 4])],
        "seed": None, "exact_maximum_arithmetic": "fractions.Fraction",
        "decimal_comparison_tolerance": 2e-15,
        "instances": [exact_squared_ratio(k, 2**16, 10000) for k in (1, 2, 8, 32, 128)]
    }, indent=2))


if __name__ == "__main__":
    main()
