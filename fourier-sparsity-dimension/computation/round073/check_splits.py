"""Verify filled spectrum and true slice sparsity across exact linear splits."""
import json
import platform
import random
from collections import Counter, defaultdict
from fractions import Fraction


def sign(value):
    return 1 - 2 * (value.bit_count() & 1)


def walsh(values):
    width = 1
    while width < len(values):
        for offset in range(0, len(values), 2 * width):
            for i in range(offset, offset + width):
                first, second = values[i], values[i + width]
                values[i], values[i + width] = first + second, first - second
        width *= 2
    return values


def basis(vectors):
    pivots = {}
    for vector in vectors:
        while vector:
            pivot = vector.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = vector
                break
            vector ^= pivots[pivot]
    return pivots


def remainder(vector, pivots):
    for pivot in sorted(pivots, reverse=True):
        if vector & (1 << pivot):
            vector ^= pivots[pivot]
    return vector


def span(vectors):
    points = [0]
    for vector in vectors:
        points += [point ^ vector for point in points]
    return points


def truth_and_support(a):
    m = 1 << a
    n = m + a - 1
    k = n - 2 * a
    t = 1 << (n + k)
    b = 2 * t
    w1_mask = (1 << n) - 1
    values = []
    for point in range(4 * t):
        x = point & w1_mask
        z = (point >> n) & ((1 << k) - 1)
        y = x & (m - 1)
        if not point & t:
            value = -1 if x == 0 else 1
        else:
            data_sign = 1 if y == 0 else sign(x & (1 << (a + y - 1)))
            value = data_sign * sign(point & b) * (-1 if z == 0 else 1)
        values.append(value)
    spectrum = {index: numerator for index, numerator in enumerate(walsh(values[:])) if numerator}
    data = [0] + [1 << (a + j) for j in range(m - 1)]
    branch_a = {x ^ selector for x in range(1 << n) for selector in (0, t)}
    branch_b = {
        b ^ selector ^ alpha ^ y ^ (z << n)
        for selector in (0, t)
        for alpha in data
        for y in range(m)
        for z in range(1 << k)
    }
    assert set(spectrum) == branch_a | branch_b
    assert len(spectrum) == 1 << (n + 2)
    assert all(spectrum[freq] == -(1 << (k + 2)) for freq in branch_a - {0, t})
    return n, k, t, b, data, values, spectrum, branch_a, branch_b


def real_slice_max(values, fixed):
    if not fixed:
        return sum(bool(c) for c in walsh(values[:]))
    representatives = {}
    for point in range(len(values)):
        label = sum(((query & point).bit_count() & 1) << i for i, query in enumerate(fixed))
        representatives.setdefault(label, point)
    kernel_basis = basis(
        point for point in range(len(values)) if all(sign(point & query) == 1 for query in fixed)
    )
    kernel = span(kernel_basis.values())
    return max(sum(bool(c) for c in walsh([values[point ^ shift] for shift in kernel]))
               for point in representatives.values())


def query_family(a, n, k, t, b, data):
    addr = [1 << j for j in range(a)]
    outer = [1 << (n + j) for j in range(k)]
    data_kernel = [
        data[j] ^ sum(data[1 << i] for i in range(a) if j & (1 << i))
        for j in range(1, len(data)) if j.bit_count() > 1
    ]
    chosen = [t] + addr + outer + data_kernel
    yield "optimal_shape", chosen, True
    yield "optimal_shape_merged", chosen + [b], True
    yield "leaving_outer_data_unqueried", [t] + addr + data_kernel + [b], True
    yield "trivial_kernel", [], True
    yield "whole_space", [1 << j for j in range(n + k + 2)], True
    rng = random.Random(7213)
    bits = n + k + 2
    for rank in range(1, bits):
        for repeat in range(5):
            vectors = []
            while len(vectors) < rank:
                candidate = rng.randrange(1, 1 << bits)
                if remainder(candidate, basis(vectors)):
                    vectors.append(candidate)
            yield f"random_rank_{rank}_{repeat}", vectors, repeat == 0


def run(a=3):
    n, k, t, b, data, values, spectrum, branch_a, branch_b = truth_and_support(a)
    m = 1 << a
    s = len(spectrum)
    rows = []
    for name, query, check_real in query_family(a, n, k, t, b, data):
        pivots = basis(query)
        frequencies_by_label = defaultdict(list)
        for freq in spectrum:
            frequencies_by_label[remainder(freq, pivots)].append(freq)
        closures = {
            label: 1 << len(basis(frequency ^ points[0] for frequency in points[1:]))
            for label, points in frequencies_by_label.items()
        }
        mass = sum(closures.values())
        assert mass >= s

        left = {remainder(frequency, pivots) for frequency in branch_a}
        right = {remainder(frequency, pivots) for frequency in branch_b}
        shared_labels = left & right
        pure_regular_labels = left - right - {0, remainder(t, pivots)}
        if pure_regular_labels:
            common_normal_size = (1 << (n + 1)) // len(left)
            for label in pure_regular_labels:
                assert closures[label] == common_normal_size
                assert len(frequencies_by_label[label]) == common_normal_size
                assert len({spectrum[freq] for freq in frequencies_by_label[label]}) == 1

        core_and_address = basis([*query, *(1 << j for j in range(a))])
        clusters = Counter(remainder(atom, core_and_address) for atom in data)
        lower_fill = (1 << k) * m * sum(1 << (size - 1) for size in clusters.values())
        assert mass >= lower_fill

        tau = k - len(basis(remainder(1 << (n + j), pivots) for j in range(k)))
        merged_lower = len(shared_labels) * ((1 << (n + 1)) // len(left)) * (1 << tau)
        assert mass >= merged_lower
        volume_overlap = sum(Fraction(1, size) for size in closures.values())
        assert volume_overlap * mass >= len(closures) ** 2
        real = real_slice_max(values, pivots.values()) if check_real else None
        if real is not None:
            assert real >= len(pure_regular_labels)
            assert real >= volume_overlap
        if name == "optimal_shape":
            assert (mass, real) == (s, m)
        rows.append({
            "split": name, "kernel_rank": len(pivots), "M": mass,
            "L": real, "L_opt_pure_label_lower": len(pure_regular_labels),
            "L_opt_volume_lower": str(volume_overlap), "merged_labels": len(shared_labels),
            "W2_kernel_rank": tau, "tail_fill_lower": lower_fill,
        })
    return {
        "python": platform.python_version(), "a": a, "m": m, "n": n, "k": k,
        "truth_entries": len(values), "s": s, "seed": 7213,
        "tested_linear_splits": len(rows),
        "tested_true_slice_maxima": sum(row["L"] is not None for row in rows),
        "named_checks": rows[:5],
        "random_checks_passed": len(rows) - 5,
    }



def ceiling_rank_variant(a=3):
    m = 1 << a
    n = m + a - 1
    k = n - 2 * a
    selector = [1 << (n + k + j) for j in range(2)]
    t_values = span(selector)
    guards = [1 << (n + k + 2 + j) for j in range(3)]
    b3 = (n + 1) // 2
    m3 = 1 << b3
    offset = n + k + 5
    addr3 = [1 << (offset + j) for j in range(b3)]
    data3 = [0] + [1 << (offset + b3 + j) for j in range(m3 - 1)]
    data = [0] + [1 << (a + j) for j in range(m - 1)]
    full_a3 = span(addr3)
    layers = [
        {x ^ code for x in range(1 << n) for code in t_values},
        {guards[0] ^ x ^ y ^ (z << n) ^ code
         for x in range(m) for y in data for z in range(1 << k) for code in t_values},
        {guards[1] ^ address ^ label ^ code
         for address in full_a3 for label in data3 for code in t_values},
        {guards[2] ^ code for code in t_values},
    ]
    spectrum = set.union(*layers)
    assert len(spectrum) == sum(map(len, layers)) == 4 * (2 ** (n + 1) + m3 ** 2 + 1)
    dim = len(basis(spectrum))
    assert dim == offset + b3 + m3 - 1

    kernel = [
        data[j] ^ sum(data[1 << i] for i in range(a) if j & (1 << i))
        for j in range(1, m) if j.bit_count() > 1
    ]
    optimal = selector + [1 << j for j in range(a)] + kernel + [
        1 << (n + j) for j in range(k)
    ] + addr3
    sampled = [optimal]
    rng = random.Random(7213)
    for requested_rank in (8, 24, 50):
        for _ in range(4):
            vectors = []
            while len(vectors) < min(requested_rank, dim):
                candidate = rng.randrange(1, 1 << dim)
                if remainder(candidate, basis(vectors)):
                    vectors.append(candidate)
            sampled.append(vectors)

    for row, query in enumerate(sampled):
        pivots = basis(query)
        by_label = defaultdict(list)
        for freq in spectrum:
            by_label[remainder(freq, pivots)].append(freq)
        mass = sum(1 << len(basis(x ^ fiber[0] for x in fiber[1:])) for fiber in by_label.values())
        space_and_a3 = basis([*query, *addr3])
        bins = Counter(remainder(label, space_and_a3) for label in data3)
        tau3 = b3 - len(basis(remainder(freq, pivots) for freq in addr3))
        assert len({remainder(freq, pivots) for freq in layers[2]}) >= len(bins) * (1 << (b3 - tau3))
        assert mass >= m3 * sum(1 << (count - 1) for count in bins.values())
        if row == 0:
            assert mass == len(spectrum)
    return {
        "a": a,
        "data_and_query_rank": dim,
        "s": len(spectrum),
        "large_branch_address_count": m3,
        "tested_sparse_linear_splits": len(sampled),
        "selected_split_M": len(spectrum),
    }


if __name__ == "__main__":
    result = run()
    result["high_rank_variant"] = ceiling_rank_variant()
    print(json.dumps(result, indent=2))
