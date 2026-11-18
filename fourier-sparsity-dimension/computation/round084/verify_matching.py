from collections import defaultdict
from fractions import Fraction
from itertools import product
from random import Random
import platform


def rank(vectors):
    pivots = {}
    for x in vectors:
        while x:
            k = x.bit_length()
            if k not in pivots:
                pivots[k] = x
                break
            x ^= pivots[k]
    return len(pivots)


def fourier_support(values):
    result = list(values)
    width = 1
    while width < len(result):
        for start in range(0, len(result), 2 * width):
            for i in range(start, start + width):
                a, b = result[i], result[i + width]
                result[i], result[i + width] = a + b, a - b
        width *= 2
    return tuple(i for i, value in enumerate(result) if value)


def matching(n, seed):
    addresses = list(range(1 << n))
    Random(seed).shuffle(addresses)
    return tuple(zip(addresses[::2], addresses[1::2]))


def explicit_support(edges, levels, policy):
    def visit(path, offset):
        if len(path) == levels:
            # All four Walsh frequencies of the two-bit bent seed are nonzero.
            return [beta << offset for beta in range(4)], offset+2, [Fraction(1)], 0
        cell = list(policy(path))
        k = len(cell).bit_length() - 1
        cursor = offset + k
        rows = []
        masses = []
        address_rank = 0
        local_directions = []
        for i, edge in enumerate(edges):
            active = [int(y in edge) for y in cell]
            if not any(active):
                continue
            local = fourier_support(active)
            local_directions.extend(local)
            guard = 1 << cursor
            data, cursor, child_masses, child_rank = visit(path + (i,), cursor+1)
            rows.extend(guard | word | (alpha << offset) for word in data for alpha in local)
            masses.extend(mu * Fraction(sum(active), len(cell)) for mu in child_masses)
            address_rank += child_rank
        return rows, cursor, masses, address_rank + rank(local_directions)

    rows, dim, masses, address_rank = visit((), 0)
    assert len(set(rows)) == len(rows)
    assert sum(masses) == 1
    assert len(rows) == 4 * sum(1 / mu for mu in masses)
    R = rank((1 << dim) | word for word in rows)
    assert R == 3 * len(masses) + address_rank
    return len(rows), R, len(masses), address_rank


def verify_recursive_spectra():
    edges = ((0, 1), (2, 4), (3, 7), (5, 6))
    checks = []
    for depth in range(1, 4):
        whole = explicit_support(edges, depth, lambda _: range(8))
        clean = explicit_support(edges, depth, lambda _: (0, 2))
        mixed = explicit_support(edges, depth, lambda _: (0, 1, 2, 3))
        # The fixed root cell has one full edge and two crossing points.
        assert whole[:2] == (4 * 4**(2 * depth), 3 * 4**depth + 3 * (4**depth - 1)//3)
        assert clean[:2] == (4 * 2**(2 * depth), 3 * 2**depth + (2**depth - 1))
        assert mixed[:2] == (4 * 10**depth, 3 * 3**depth + (3**depth - 1))
        checks.extend([whole, clean, mixed])

    def compensated(path):
        if not path:
            return (0, 1, 2, 3)
        return (0, 2) if path[0] == 0 else (0,)

    assert explicit_support(edges, 2, compensated)[0] == 4 * 4**2
    return len(checks) + 1


def codewords(basis, offset):
    values = [offset]
    for v in basis:
        values += [u ^ v for u in values]
    return values


def address_subgroup_checks():
    n, L, N = 3, 2, 4
    edges = ((0, 1), (2, 4), (3, 7), (5, 6))
    which = {y: i for i, edge in enumerate(edges) for y in edge}
    dim = n * (1 + N)
    rng = Random(84021)
    for case in range(320):
        span_dim = case % 9
        basis = []
        while len(basis) < span_dim:
            x = rng.randrange(1, 1 << dim)
            if rank(basis + [x]) == len(basis) + 1:
                basis.append(x)
        fibres = codewords(basis, rng.randrange(1 << dim))
        addresses = {0: {x & 7 for x in fibres}}
        for x in fibres:
            child = 1 + which[x & 7]
            addresses.setdefault(child, set()).add((x >> (n * child)) & 7)
        charged = 0
        previous_mask = 0
        last_rank = 0
        for node, possible in sorted(addresses.items()):
            assert len(possible) and len(possible) & (len(possible) - 1) == 0
            local_codim = n - (len(possible).bit_length() - 1)
            previous_mask |= 7 << (n * node)
            new_rank = rank(x & previous_mask for x in basis)
            available = n + last_rank - new_rank
            assert 0 <= local_codim <= available
            charged += local_codim
            last_rank = new_rank
        assert charged <= dim - span_dim

    # Exact selected-leaf certificate: root on edge (0,1), chosen child's
    # address on edge (2,4). Their affine directions are 1 and 6.
    free = [1, 6 << 3] + [1 << j for j in range(6, dim)]
    residue = codewords(free, 2 << 3)
    assert { (which[x & 7], which[(x >> 3) & 7]) for x in residue } == {(0, 1)}
    assert dim - rank(free) == (n - 1) * L

    # The linear span of all leaf indicator Fourier supports is the annihilator
    # of the period space of the entire leaf map.
    band = [fourier_support(int(x in edge) for x in range(8)) for edge in edges]
    support = {a | (b << (n * (i+1))) for i in range(N) for j in range(N)
               for a in band[i] for b in band[j]}
    assert rank(support) == dim
    return 320, dim


def shared_address_checks():
    edges = ((0, 1), (2, 4), (3, 7), (5, 6))
    bands = [fourier_support(int(y in edge) for y in range(8)) for edge in edges]
    checks = []
    for L in range(1, 4):
        address_bits = 3 * L
        frequencies = []
        for leaf, path in enumerate(product(range(4), repeat=L)):
            offset = address_bits + 3 * leaf
            data = [(1 | (freq << 1)) << offset for freq in range(4)]
            for address in product(*(bands[i] for i in path)):
                alpha = sum(freq << (3 * j) for j, freq in enumerate(address))
                frequencies.extend(alpha | word for word in data)
        assert len(frequencies) == 4 * 4**(2*L)
        R = rank(frequencies)
        assert R == 3 * 4**L + address_bits
        checks.append((len(frequencies), R))
    return checks


def compensated_larger_certificate():
    n = 10
    N = 1 << (n-1)
    b = n // 2
    M = 1 << b
    edges = matching(n, 84022)
    which = {y: i for i, edge in enumerate(edges) for y in edge}
    for high in range(1 << (n-b)):
        cell = tuple((high << b) | low for low in range(M))
        if len({which[x] for x in cell}) == M:
            break
    else:
        raise AssertionError("no empty matched cell")
    root = set(range(N))
    f = sum(x in root and y in root for x, y in edges)
    t = N - 2*f
    assert f and t
    B = f*M + t*(M//2)
    assert B == N*M//2
    masses = {Fraction(2, N*M), Fraction(1, N*(M//2))}
    assert masses == {Fraction(1, B)}
    parent_D = 3*N*N + n*(N+1) + 1
    child_D = 3*B+(n-1) + f*b + t*(b-1) + 1
    assert child_D * N*N > parent_D * B
    return n, high, f, t, B, (parent_D, 4*N**4), (child_D, 4*B**2)


if __name__ == "__main__":
    print("Python", platform.python_version())
    print("explicit recursive guarded spectra:", verify_recursive_spectra())
    print("mixed affine fibres and dimension of all leaf normal span:", address_subgroup_checks())
    print("compensated two-generation certificate:", compensated_larger_certificate())
    print("shared addresses per depth, (s,R):", shared_address_checks())
