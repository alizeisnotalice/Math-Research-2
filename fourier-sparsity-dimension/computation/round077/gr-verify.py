from collections import Counter
from itertools import combinations
import platform


def fwht(values):
    values = list(values)
    step = 1
    while step < len(values):
        for base in range(0, len(values), 2 * step):
            for index in range(base, base + step):
                left, right = values[index], values[index + step]
                values[index], values[index + step] = left + right, left - right
        step *= 2
    return values


def check_spectrum(values, k, address_bits, expected):
    n = len(values).bit_length() - 1
    scale = 1 << (n - k)
    blocks = {}
    for frequency, value in enumerate(fwht(values)):
        assert value % scale == 0
        coefficient = value // scale
        if not coefficient:
            continue
        beta = frequency >> address_bits
        if beta == 0:
            assert coefficient & 1
        else:
            assert not coefficient & 1
        blocks.setdefault(beta, []).append(coefficient)
    assert len(blocks[0]) == 1 << address_bits
    assert set(blocks) == {0, *expected}
    for beta, (length, magnitude) in expected.items():
        assert len(blocks[beta]) == length
        assert all(abs(value) == magnitude for value in blocks[beta])


def check_single_layer(k, K, rank):
    address_bits = k + 1
    u_mask = (1 << k) - 1
    tail_selector_mask = (1 << (K - rank)) - 1
    values = []
    for index in range(1 << (address_bits + rank)):
        u = index & u_mask
        a = (index >> k) & 1
        z = index >> address_bits
        negative = (u == 0) if a == 0 else (z == 0 and u & tail_selector_mask == 0)
        values.append(1 - 2 * negative)
    expected = {
        beta: (1 << (K - rank + 1), 1 << (k - K))
        for beta in range(1, 1 << rank)
    }
    check_spectrum(values, k, address_bits, expected)


def check_three_layer_boolean_realization():
    k = 6
    depths = (5, 4, 3)
    address_bits = k + 1
    addresses = [None] * (1 << address_bits)
    cursor = 0
    beta_index = 0
    expected = {}
    for K in reversed(depths):
        width = 1 << (address_bits - K)
        assert cursor % width == 0
        for _ in range(K):
            addresses[cursor : cursor + width] = [beta_index] * width
            expected[1 << beta_index] = (1 << K, 1 << (k - K))
            cursor += width
            beta_index += 1
    assert cursor == 100
    values = []
    for z in range(1 << beta_index):
        for address, selected in enumerate(addresses):
            negative = (address == cursor) if selected is None else ((z >> selected) & 1)
            values.append(1 - 2 * negative)
    check_spectrum(values, k, address_bits, expected)
    return sum(length for length, _ in expected.values())


def check_anchor_uniqueness():
    tested = 0
    layer1 = list(range(1, 8))
    for b2 in range(8):
        layer2 = [b2 ^ (u << 3) for u in range(1, 4)]
        for b3 in range(32):
            layer3 = [b3 ^ (1 << 5)]
            multiplicities = Counter(
                u ^ v for u, v in combinations([0, *layer1, *layer2, *layer3], 2)
            )
            for earlier, later, anchor in (
                (layer1, layer2, b2), (layer1, layer3, b3), (layer2, layer3, b3)
            ):
                assert all(
                    multiplicities[u ^ v] == 1 for u in earlier for v in later if u != anchor
                )
            assert any(u not in (b2, b3) for u in layer1)
            assert any(u != b3 for u in layer2)
            tested += 1
    return tested


def main():
    print("Python", platform.python_version())
    single_count = 0
    for k in range(2, 7):
        for K in range(1, k):
            for rank in range(1, K + 1):
                check_single_layer(k, K, rank)
                single_count += 1
    print("single-layer Boolean equality instances:", single_count)
    print("three-layer Boolean construction, actual q:", check_three_layer_boolean_realization())
    print("three-layer anchor patterns checked:", check_anchor_uniqueness())


if __name__ == "__main__":
    main()
