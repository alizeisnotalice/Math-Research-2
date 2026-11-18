"""Exact weighted spectral-profile and gluing tests for affine terminal-character rows."""

from fractions import Fraction
import json
import platform
import random


def walsh(values):
    result = list(values)
    width = 1
    while width < len(result):
        for start in range(0, len(result), 2 * width):
            for index in range(start, start + width):
                a, b = result[index], result[index + width]
                result[index], result[index + width] = a + b, a - b
        width *= 2
    return {index: Fraction(coefficient, len(result))
            for index, coefficient in enumerate(result) if coefficient}


def rank(rows):
    pivots = {}
    for row in rows:
        while row:
            index = row.bit_length() - 1
            if index not in pivots:
                pivots[index] = row
                break
            row ^= pivots[index]
    return len(pivots)


def bit(word, index):
    return (word >> index) & 1


def label_function(blocks):
    points = sum(map(len, blocks))
    selector_dimension = points.bit_length() - 1
    assert points == 1 << selector_dimension
    selected_labels = [next(i for i, block in enumerate(blocks) if x in block) for x in range(points)]
    samples = [1 - 2 * bit(y, selected_labels[x])
               for y in range(1 << len(blocks)) for x in range(points)]
    return walsh(samples), len(samples)


def audit_density(blocks, records):
    spectrum, points = label_function(blocks)
    full_rank = rank(frequency | points for frequency in spectrum)
    normal_rank = full_rank - len(blocks)
    delta = len(blocks) - normal_rank
    assert delta >= 1 and sum(map(abs, spectrum.values())) == len(blocks)
    entries = list(spectrum.items())
    tested = 0
    for mask in records(len(entries)):
        subset = [index for pos, (index, _) in enumerate(entries) if bit(mask, pos)]
        r = rank(index | points for index in subset)
        if r == full_rank:
            continue
        mass = sum(abs(spectrum[index]) for index in subset)
        assert 2 * mass <= r + delta - 1
        assert full_rank + 1 - r <= 2 * (len(blocks) - mass)
        tested += 1
    return {"size": len(spectrum), "full_augmented_rank": full_rank, "delta": delta,
            "non_full_rank_subsets_checked": tested}


def nonparallel_multilift(four_lines):
    chosen = [point for point, _ in four_lines]
    labels = [0] * 8
    for index, (first, second) in enumerate(four_lines, 1):
        labels[second] = index
    values = [1 - 2 * bit(y, labels[x]) for y in range(32) for x in range(8)]
    spectrum = walsh(values)
    local_shared = walsh([1, 0])
    global_D = rank(index | len(values) for index in spectrum) + 1
    global_L1 = sum(map(abs, spectrum.values()))
    assert len(local_shared) == 2 and set(chosen) == {0, 2, 3, 5}
    assert global_D <= 2 * global_L1
    return {"s": len(spectrum), "D": global_D, "L1": str(global_L1), "local_shared_lifts": 2}


def nonreciprocal_coefficient():
    values = []
    for word in range(128):
        w, z1, z2 = bit(word, 0), bit(word, 1), bit(word, 2)
        if not w:
            label = 1 if z1 else 0
        elif not z1:
            label = 2
        else:
            label = 3 if z2 else 0
        values.append(1 - 2 * bit(word, 3 + label))
    actual = walsh(values)[8]
    assert actual == Fraction(3, 8)
    return str(actual)


def recurrence_counterexample():
    data = {}
    for dimension in (2, 3, 4):
        values = [1] * ((1 << dimension) - 1) + [-1]
        spectrum = walsh(values)
        D = rank(index | len(values) for index in spectrum) + 1
        L1 = sum(map(abs, spectrum.values()))
        residual_L1 = L1 - abs(spectrum[0])
        data[dimension] = {"D": D, "s": len(spectrum), "L1": str(L1),
                           "outside_DC_L1": str(residual_L1)}
    assert data[3]["D"] == 2 * Fraction(data[3]["L1"])
    assert data[3]["D"] - 1 > 2 * Fraction(data[3]["outside_DC_L1"])
    assert data[4]["D"] > 2 * Fraction(data[4]["L1"])
    return data


def main():
    four_lines = [(0, 1), (2, 4), (3, 7), (5, 6)]
    unequal = [(0,), (1,)] + four_lines[1:]
    seed = 510426
    rng = random.Random(seed)
    result = {
        "python": platform.python_version(),
        "seed": seed,
        "equal_mass_nonparallel_exhaustive": audit_density(four_lines, lambda n: range(1 << n)),
        "unequal_mass_nonparallel_sampled": audit_density(
            unequal, lambda n: [rng.randrange(1 << n) for _ in range(16384)]
        ),
        "nonparallel_outer_multi_lift": nonparallel_multilift(four_lines),
        "nonreciprocal_global_amplitude_after_valid_step": nonreciprocal_coefficient(),
        "recursion_counterexample": recurrence_counterexample(),
    }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
