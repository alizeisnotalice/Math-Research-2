"""Verify data-frequency shearing and the nonlinear-offset obstruction."""

import json
import platform


def walsh(values):
    values = list(values)
    width = 1
    while width < len(values):
        for start in range(0, len(values), 2 * width):
            for index in range(start, start + width):
                a, b = values[index], values[index + width]
                values[index], values[index + width] = a + b, a - b
        width *= 2
    return values


def bit(word, index):
    return (word >> index) & 1


def character(bit_value):
    return 1 - 2 * bit_value


def nonlinear_case(m):
    """Root x0+OR(w), two occurrences of x1, left target only at w=0."""
    def value(private, w):
        bias = w != 0
        root = bit(private, 0) ^ bias
        terminal = character(bit(private, 1))
        if not root:
            return terminal if w == 0 else 1
        return 1 if w == 0 else terminal

    private_size, shared_size = 4, 1 << m
    truth = [value(private, w) for w in range(shared_size) for private in range(private_size)]
    assert all(value(x, w) == (character(bit(x, 1)) if not bit(x, 0) else 1)
               for x in range(private_size) for w in range(shared_size))
    actual = [index for index, coefficient in enumerate(walsh(truth)) if coefficient]
    incidence = sum(index & 3 in (2, 3) for index in actual)
    constrained = []
    for w in range(shared_size):
        private = int(w != 0)
        constrained.append((value(private, w) - value(private ^ 2, w)) // 2)
    conditioned_support = sum(c != 0 for c in walsh(constrained))
    assert constrained == [1] + [0] * (shared_size - 1)
    assert actual == [0, 1, 2, 3] and incidence == 2 and conditioned_support == shared_size
    return {"data_dimension": m, "true_sparsity": len(actual), "true_band_incidence": incidence,
            "isolated_signal_sparsity": conditioned_support, "unpaid_claim": shared_size}


def affine_case(m):
    """The target x1 is isolated using one shared cut and one private shear."""
    def value(private, w):
        terminal = character(bit(private, 1))
        if bit(w, 0):
            return terminal
        root = bit(private, 0) ^ bit(w, 1)
        if root:
            return terminal
        return terminal if not w else 1

    conditional_size = 1 << (m - 1)
    values = [value(x, w) for w in range(1 << m) for x in range(4)]
    original_spectrum = walsh(values)
    original_incidence = sum(coefficient != 0 and index & 3 in (2, 3)
                             for index, coefficient in enumerate(original_spectrum))
    conditional_signal = [
        (value(bit(w, 1), w) - value(bit(w, 1) ^ 2, w)) // 2
        for w in range(0, 1 << m, 2)
    ]
    direct = walsh(conditional_signal)
    # On w0=0 and x0=w1, the private x0 frequency becomes a translation in w1.
    shear = [0] * conditional_size
    for frequency, coefficient in enumerate(original_spectrum):
        if frequency & 3 not in (2, 3):
            continue
        beta = frequency >> 2
        beta_restricted = beta >> 1
        sheared_frequency = beta_restricted ^ bit(frequency, 0)
        shear[sheared_frequency] += coefficient
    assert shear == [8 * value for value in direct]
    assert conditional_signal == [1] + [0] * (conditional_size - 1)
    assert original_incidence >= conditional_size
    return {"data_dimension": m, "true_band_incidence": original_incidence,
            "isolated_signal_sparsity": sum(c != 0 for c in direct),
            "data_frequencies_shear_one_for_one": True}


def main():
    result = {"python": platform.python_version(),
              "nonlinear_failures": [nonlinear_case(m) for m in range(2, 7)],
              "affine_shear_checks": [affine_case(m) for m in range(2, 7)]}
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
