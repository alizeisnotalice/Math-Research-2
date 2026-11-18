"""Exact arithmetic checks of the graph realization and the all-quotient arguments."""
from fractions import Fraction
from itertools import product
import json
import platform


def walsh(values):
    w = list(values)
    k = 1
    while k < len(w):
        for off in range(0, len(w), 2*k):
            for i in range(off, off+k):
                a, b = w[i], w[i+k]
                w[i], w[i+k] = a + b, a - b
        k *= 2
    return w


def graph(n, a, b):
    u_max = 1 << n
    size = 2*u_max**2
    F = []
    for i in range(size):
        t = i & 1
        x = (i >> 1) & (u_max-1)
        y = i >> (n+1)
        bit = ((x & y).bit_count() & 1) ^ (((b >> x) if t else (a >> y)) & 1)
        F.append(1 - 2*bit)
    expected = set()
    w = walsh(F)
    for u, v in product(range(u_max), repeat=2):
        f_bit = ((a >> u) ^ (b >> v)) & 1
        xi = (v << (n+1)) | (u << 1) | f_bit
        signed_value = u_max*2 * (-1)**(((u & v).bit_count() + ((a >> u) & 1)) % 2)
        assert w[xi] == signed_value
        expected.add(xi)
    assert expected == {i for i, z in enumerate(w) if z}
    return expected


def linear_subspaces(n):
    spaces = {frozenset([0])}
    frontier = set(spaces)
    while frontier:
        new = set()
        for H in frontier:
            for x in range(1 << n):
                if x not in H:
                    new.add(H | {z ^ x for z in H})
        frontier = new - spaces
        spaces.update(frontier)
    return spaces


def affine_hull(A):
    anchor = next(iter(A))
    dirs = {0}
    for u in A:
        v = u ^ anchor
        if v not in dirs:
            dirs |= {z ^ v for z in dirs}
    return {anchor ^ z for z in dirs}


def cube_count(A, k):
    # Ordered affine maps with an arbitrary basepoint and k arbitrary directions.
    if not k:
        return len(A)
    if not A:
        return 0
    ambient = 1 << max(A).bit_length()
    return sum(cube_count(A & {x ^ v for x in A}, k-1) for v in range(ambient))


def all_quotient_check(A, n):
    total = 0
    holes = 0
    full_cubes = cube_count(A, 3)
    for H in linear_subspaces(n):
        bands = {}
        for x in A:
            key = min(x ^ h for h in H)
            bands.setdefault(key, set()).add(x)
        sizes = [len(affine_hull(B)) for B in bands.values()]
        density_sum = sum((Fraction(1, z) for z in sizes), Fraction())
        # The overlap lower bound max(1, average load) suffices for the proof.
        E = max(Fraction(1), density_sum) * sum(sizes) / len(A)
        assert len(bands)**2 <= E*len(A)
        boxes = [cube_count(B, 3) for B in bands.values()]
        for B, hull, cnt in zip(bands.values(), sizes, boxes):
            assert cnt * hull**4 >= len(B)**8
        assert full_cubes >= sum(boxes)
        # Square both sides of (1) to clear the half exponents for k=3.
        assert full_cubes**2 * 2**16 * E**11 >= len(A)**5
        total += 1
        holes += bool(sum(sizes) > len(A))
    return {"all_kernels": total, "with_fibre_holes": holes, "labelled_3_cubes": full_cubes}


def rational_certificates():
    n = 24
    ambient = 2*n + 1
    log_candidates = (2*n+1)*(n+1)
    samples = 1 << (n//2)
    assert 5**3 < 2**7
    assert 3 * log_candidates < samples
    assert 6**4 > 2 * 5**4
    k = 2500
    gap = (2*n)*(1 + k//2) - (ambient*(ambient+1) + (n-1)*(k+1))
    assert gap == 75
    # Exact Chernoff bound cubed, including the union over affine functions,
    # avoiding a nonintegral power in the result.
    union_bound_cubed = Fraction(2**(3*log_candidates), 1) * Fraction(125, 256)**samples
    assert union_bound_cubed < 1
    assert 3*log_candidates - samples == -421
    return {"n": n, "log2_count_plus_targets_upper": log_candidates,
            "independent_fibres_lower": samples,
            "log2_failure_upper_after_cubing": 3*log_candidates-samples,
            "cube_order": k, "log2_cube_gap_lower": gap,
            "direct_exponent_numerator": "2*m-210", "direct_exponent_denominator": 209}


def main():
    checks = {}
    for n in (1, 2):
        count = 0
        for a, b in product(range(1 << (1 << n)), repeat=2):
            graph(n, a, b)
            count += 1
        checks[n] = count
    A = graph(2, 1, 8)
    tensor = {x | (y << 5) for x, y in product(A, repeat=2)}
    assert cube_count(tensor, 2) == cube_count(A, 2)**2
    print(json.dumps({"python": platform.python_version(),
                      "all_pairs_of_small_truth_tables_checked": checks,
                      "five_variable_graph_all_frequency_kernels": all_quotient_check(A, 5),
                      "tensor_labelled_2_cubes": cube_count(tensor, 2),
                      "symbolic_large_seed_and_exponent": rational_certificates()}, indent=2))


if __name__ == "__main__":
    main()
