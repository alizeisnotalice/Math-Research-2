"""Check all Fourier bands after one inactive address block is split into singletons."""

from fractions import Fraction
from itertools import product
import json
from math import isqrt
import platform


def fwht(values):
    values = list(values)
    width = 1
    while width < len(values):
        for start in range(0, len(values), 2 * width):
            for i in range(start, start + width):
                a, b = values[i], values[i + width]
                values[i], values[i + width] = a + b, a - b
        width *= 2
    return [Fraction(value, len(values)) for value in values]


def gf2_rank(vectors):
    pivots = {}
    for vector in vectors:
        while vector:
            pivot = vector.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = vector
                break
            vector ^= pivots[pivot]
    return len(pivots)


def exact_sqrt(value):
    n, d = isqrt(value.numerator), isqrt(value.denominator)
    assert n * n == value.numerator and d * d == value.denominator
    return Fraction(n, d)


def parity(x):
    return x.bit_count() % 2


def make_cells(address_signs):
    count = len(address_signs)
    split_address = next(i for i in range(count) if not address_signs[i])
    cells = []
    label = 1  # label zero denotes the shared data basis vector
    for address, sign in enumerate(address_signs):
        if sign:
            outputs = []
            for point in range(8):
                data_label = label if point & 4 else 0
                bent_sign = 1 - 2 * ((point & 3) == 3)
                outputs.append((point, data_label, sign * bent_sign))
            cells.append((address, outputs))
            label += 1
        elif address == split_address:
            for point in range(8):
                cells.append((address, [(point, label, 1)]))
                label += 1
        else:
            cells.append((address, [(point, label, 1) for point in range(8)]))
            label += 1
    assert label == count + 8
    return cells, split_address, label


def check(address_signs, use_signed):
    m = len(address_signs)
    r = m.bit_length() - 1
    u = sum(sign != 0 for sign in address_signs)
    assert m == 1 << r
    selector_size = 8 * m
    selector_dim = r + 3
    cells, split_address, number_labels = make_cells(address_signs)
    assert len(cells) == m + 7

    direct_signals = [[0] * selector_size for _ in range(number_labels)]
    lifted = [[Fraction(0)] * selector_size for _ in range(number_labels)]
    physical_partition = set()
    shared_phases = [0] * selector_size
    common_tangent = set(range(selector_size))
    local_energy = []
    local_budgets = []
    mean_shared_abs_mass = []
    sum_local_D = 0
    lifted_local_beta_count = []
    for address, outputs in cells:
        points = {(address << 3) | local for local, _, _ in outputs}
        assert physical_partition.isdisjoint(points)
        physical_partition.update(points)
        reference = next(iter(points))
        common_tangent.intersection_update(reference ^ point for point in points)
        p_i = Fraction(len(outputs), selector_size)
        labels = sorted({label for _, label, _ in outputs})
        partials = {
            label: fwht([sign if entry_label == label else 0 for _, entry_label, sign in outputs])
            for label in labels
        }
        local_dimension = len(outputs).bit_length() - 1
        local_augmented = [
            direction | (1 << (local_dimension + label)) | (1 << (local_dimension + number_labels))
            for label, coefficients in partials.items()
            for direction, coefficient in enumerate(coefficients) if coefficient
        ]
        s_i = len(local_augmented)
        D_i = gf2_rank(local_augmented) + 1
        sum_local_D += D_i
        assert (D_i, s_i) == ((6, 16) if 0 in labels else (2, 1))

        for local_point, label, sign in outputs:
            direct_signals[label][(address << 3) | local_point] = sign

        beta_spectrum = partials.get(0, [])
        if beta_spectrum:
            assert set(map(abs, beta_spectrum)) == {Fraction(1, 4)}
            beta_count = sum(value != 0 for value in beta_spectrum)
            lifted_local_beta_count.append(beta_count)
            energy = sum(value**2 for value in beta_spectrum)
            local_energy.append(p_i * energy)
            local_budgets.append(p_i * Fraction(D_i, s_i)**2 * beta_count)
            mean_shared_abs_mass.append(sum(map(abs, beta_spectrum)))

        for selector_frequency in range(selector_size):
            address_frequency = selector_frequency >> 3
            local_frequency = selector_frequency & 7
            address_phase = 1 - 2 * parity(address & address_frequency)
            for label, spectrum in partials.items():
                if len(outputs) == 1:
                    private_phase = 1 - 2 * parity(local_frequency & outputs[0][0])
                    lifted[label][selector_frequency] += p_i * spectrum[0] * address_phase * private_phase
                else:
                    contribution = p_i * spectrum[local_frequency] * address_phase
                    lifted[label][selector_frequency] += contribution
                    if label == 0:
                        shared_phases[selector_frequency] += 1 if contribution > 0 else -1

    assert physical_partition == set(range(selector_size))
    assert common_tangent == {0}, "the ambient shared tangent still contains selector directions"
    assert lifted_local_beta_count == [8] * u
    direct = list(map(fwht, direct_signals))
    assert direct == lifted
    energy_beta = sum(local_energy)
    assert energy_beta == Fraction(u, 2 * m)
    assert sum(value**2 for value in direct[0]) == energy_beta
    shared_s = sum(value != 0 for value in direct[0])
    expected_shared = 64 if use_signed else 32
    assert shared_s == expected_shared
    expected_private = m * m + 7 * m * u + 63 * m
    actual_private = sum(value != 0 for band in direct[1:] for value in band)
    assert actual_private == expected_private

    entire_support = [
        selector_frequency | (1 << (selector_dim + label))
        for label, band in enumerate(direct)
        for selector_frequency, coefficient in enumerate(band) if coefficient
    ]
    whole_rank = gf2_rank(entire_support)
    whole_D = gf2_rank(
        field_frequency | (1 << (selector_dim + number_labels))
        for field_frequency in entire_support
    ) + 1
    assert len(entire_support) == expected_private + expected_shared
    assert whole_rank == m + r + 11 and whole_D == m + r + 12
    assert sum_local_D == 4 * u + 2 * m + 14
    projected_rank = gf2_rank(
        (1 << label) | (1 << number_labels)
        for label, band in enumerate(direct)
        if any(band)
    )
    normal_rank = selector_dim
    B_rank = sum_local_D - (normal_rank + projected_rank + 1)
    assert B_rank == 4 * u + m + 2 - r

    max_phase_sum = max(map(abs, shared_phases))
    lambda_energy = Fraction(8 * m, u)
    lambda_signed = sum(mean_shared_abs_mass)**2 / (energy_beta * max_phase_sum**2)
    lambda_use = lambda_signed if use_signed else lambda_energy
    assert lambda_use <= shared_s
    root_term = exact_sqrt(lambda_use * sum(local_budgets))
    rhs = 3 * u - root_term
    assert B_rank >= rhs
    if use_signed:
        assert lambda_use == Fraction(320, 9) and max_phase_sum == 3 and rhs == 10
    else:
        assert lambda_use == Fraction(32, 3) and max_phase_sum == u and rhs == 3 * u - 3

    return {
        "address_dimension": r,
        "address_count": m,
        "participant_count": u,
        "split_address": split_address,
        "outer_blocks": len(cells),
        "data_labels": number_labels,
        "selector_tangent_dimension_after_split": 0,
        "outer_normal_rank": normal_rank,
        "exact_compared_fourier_coefficients": selector_size * number_labels,
        "shared_sparsity": shared_s,
        "private_sparsity": actual_private,
        "s": len(entire_support),
        "d": whole_rank,
        "D": whole_D,
        "rank_budget_lhs": B_rank,
        "rank_budget_rhs": str(rhs),
        "lambda_energy": str(lambda_energy),
        "lambda_signed": str(lambda_signed),
        "max_signed_sum": max_phase_sum,
    }


def main():
    results = {
        "python": platform.python_version(),
        "equal_energy": [check([1 if (w & 3) != 3 else 0 for w in range(1 << r)], False)
                         for r in (2, 3, 5)],
        "signed": check([1, 1, 1, -1, 1, 0, 0, 0], True),
    }
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
