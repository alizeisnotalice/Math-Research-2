from fractions import Fraction
from itertools import product
import platform


def walsh(values):
    result = list(values)
    width = 1
    while width < len(result):
        for start in range(0, len(result), 2 * width):
            for index in range(start, start + width):
                left, right = result[index], result[index + width]
                result[index] = left + right
                result[index + width] = left - right
        width *= 2
    return result


def odd_ternary_profile_checks():
    n = 4
    N = 1 << n
    shift = N - 1
    block_masks = (
        (1, 2, 4, 8), (3, 4, 8), (3, 12), (5, 10), (9, 6), (7, 8)
    )
    strict_checks = 0
    boundary_crossing = 0
    for initial in product((-1, 0, 1), repeat=N // 2):
        values = list(initial) + [-value for value in reversed(initial)]
        assert all(values[x ^ shift] == -values[x] for x in range(N))
        coeff = walsh(values)
        assert all(value == 0 for freq, value in enumerate(coeff) if freq.bit_count() % 2 == 0)
        for blockset in block_masks:
            blocks = [
                walsh([value if freq and freq & mask == freq else 0
                       for freq, value in enumerate(coeff)])
                for mask in blockset
            ]
            bridge = [N * values[x] - sum(block[x] for block in blocks) for x in range(N)]
            infinity = max(map(abs, bridge))
            if 4 * infinity < N:
                strict_checks += 1
                assert all(bridge_value == 0 for bridge_value in bridge)
                active = [block for block in blocks if any(block)]
                if len(active) == 2:
                    assert all(abs(value) == N // 2 for block in active for value in block)
                else:
                    assert len(active) <= 1
            elif 4 * infinity == N and any(bridge):
                boundary_crossing += 1
    return N, len(block_masks), strict_checks, boundary_crossing


def modulated_bridge_examples():
    outcomes = []
    for m in range(9):
        delta = Fraction(1, 1 << m)
        count = 1 << (m + 3)
        ternary = []
        for index in range(count):
            x = 1 - 2 * (index & 1)
            y = 1 - 2 * ((index >> 1) & 1)
            z = 1 - 2 * ((index >> 2) & 1)
            value = (x + y + z + x * y * z) // 4 if index >> 3 == 0 else 0
            ternary.append(value)
        F = [1 - h * h + h for h in ternary]
        assert all(value in (-1, 1) for value in F)
        assert all((value - F[index ^ 7]) // 2 == ternary[index]
                   for index, value in enumerate(F))
        num = walsh(ternary)
        bridge_nums = [a for beta, a in enumerate(num) if beta & 7 == 7]
        assert len(bridge_nums) == 1 << m
        assert all(abs(a) == 2 for a in bridge_nums)
        assert sum(Fraction(a, count) ** 2 for a in bridge_nums) == delta / 16
        outcomes.append(str(delta))
    return outcomes


def rare_two_block_checks():
    observed = []
    for exponent in range(1, 13):
        p = Fraction(1, 1 << exponent)
        outcomes = ((-1, p / 2), (0, 1 - p), (1, p / 2))
        second_moment = 0
        for a, p_a in outcomes:
            for b, p_b in outcomes:
                H = max(-1, min(1, a + b))
                second_moment += p_a * p_b * (H - (1 - p / 2) * (a + b)) ** 2
            conditional_average = sum(q * max(-1, min(1, a + b)) for b, q in outcomes)
            assert conditional_average == (1 - p / 2) * a
        assert second_moment == p * p * (1 - p) / 2
        assert p < 2 * (1 - p)
        observed.append(str(second_moment))
    return observed


if __name__ == "__main__":
    print("Python", platform.python_version())
    print("four-variable odd ternary truth tables:", 3**8)
    print("N, partitions, strict checks, crossing boundary checks:", odd_ternary_profile_checks())
    print("gated bridge equality densities:", modulated_bridge_examples())
    print("two rare blocks, exact squared errors:", rare_two_block_checks())
