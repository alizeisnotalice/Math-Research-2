from functools import lru_cache
from random import Random
import hashlib
import json
import platform


def matching(n, seed):
    addresses = list(range(1 << n))
    Random(seed).shuffle(addresses)
    return list(zip(addresses[::2], addresses[1::2]))


def rank(vectors):
    pivots = {}
    for x in vectors:
        while x:
            key = x.bit_length()
            if key not in pivots:
                pivots[key] = x
                break
            x ^= pivots[key]
    return len(pivots)


@lru_cache(None)
def spectral_rows(values):
    values = list(values)
    width = 1
    while width < len(values):
        for start in range(0, len(values), 2 * width):
            for i in range(start, start + width):
                a, b = values[i], values[i + width]
                values[i], values[i + width] = a + b, a - b
        width *= 2
    return tuple(i for i, value in enumerate(values) if value)


def query_echelon(equations):
    pivots = {}
    for d, phase, answer in equations:
        while d:
            key = d.bit_length()
            if key not in pivots:
                pivots[key] = (d, phase, answer)
                break
            r, q, s = pivots[key]
            d, phase, answer = d ^ r, phase ^ q, answer ^ s
        else:
            raise AssertionError("dependent data equations")
    return tuple(sorted(pivots.items(), reverse=True))


def reduce_label(label, echelon):
    phase = sign = 0
    for key, (other, q, s) in echelon:
        if label >> (key - 1) & 1:
            label ^= other
            phase ^= q
            sign ^= s
    return label, phase, sign


def direct_small_check(edges, cell, equations):
    N = len(edges)
    M = len(cell)
    data_bits = 3 * N
    eqs = query_echelon(equations)
    index = {point: rank for rank, point in enumerate(cell)}
    F = {}
    occurrence_ranks = 0
    occurrence_count = 0
    B = f = t = 0
    for i, original_edge in enumerate(edges):
        points = [index[point] for point in original_edge if point in index]
        if not points:
            continue
        B += 1
        f += len(points) == 2
        t += len(points) == 1
        local = {}
        for local_freq in range(4):
            d = (1 | (local_freq << 1)) << (3 * i)
            beta, phase, sign = reduce_label(d, eqs)
            values = local.setdefault(beta, [0] * M)
            for u in points:
                negative = ((phase & u).bit_count() & 1) ^ sign ^ (local_freq == 3)
                values[u] += 1 - 2 * negative
        active = [beta for beta, vals in local.items() if any(vals)]
        for beta in active:
            dest = F.setdefault(beta, [0] * M)
            for u in points:
                assert dest[u] == 0
                dest[u] += local[beta][u]
        occurrence_ranks += rank(active)
        occurrence_count += len(active)
    r = rank(F)
    s = 0
    frequencies = []
    for beta, values in F.items():
        band = spectral_rows(tuple(values))
        assert band
        for alpha in band:
            frequencies.append(beta | (alpha << data_bits))
        s += len(band)
    D = rank(frequencies) + 1
    L = 3 * B - r
    J = occurrence_ranks - r
    dup = occurrence_count - len(F)
    w = len(equations)
    m = M.bit_length() - 1
    assert D <= 3 * B + m + 1 - L
    assert 0 <= L <= 2 * w
    assert 3 * dup <= 4 * J
    if M == 2 * N:
        assert w <= L
        assert s >= 4 * N * N - 2 * N * L
    else:
        assert 3 * s >= 6 * M * (f + 2 * t) - 8 * M * L
    return D, s


def verify_small_mixed():
    n = 3
    N = 1 << (n - 1)
    edges = matching(n, 8301)
    assert len({x ^ y for x, y in edges}) > 1
    guard = sum(1 << (3 * i) for i in range(N))
    masks = [d for d in range(1, 1 << (3 * N)) if (d & guard).bit_count() % 2 == 0]
    count = 0
    for d in masks:
        for slope in range(1 << n):
            for bit in (0, 1):
                direct_small_check(edges, range(1 << n), [(d, slope, bit)])
                count += 1
    rng = Random(8302)
    sampled = [0, 0, 0]
    for trial in range(2500):
        a = trial % 3
        normals = rng.sample(range(1, 1 << n), a)
        while rank(normals) < a:
            normals = rng.sample(range(1, 1 << n), a)
        offset = rng.randrange(1 << n)
        equations = []
        w = trial % 4
        cell = [y for y in range(1 << n)
                if all((normal & (y ^ offset)).bit_count() % 2 == 0 for normal in normals)]
        m = len(cell).bit_length() - 1
        data_masks = rng.sample(masks, w)
        while rank(data_masks) < w:
            data_masks = rng.sample(masks, w)
        for d in data_masks:
            equations.append((d, rng.randrange(1 << m), rng.randrange(2)))
        direct_small_check(edges, cell, equations)
        sampled[a] += 1
    return count, sampled


def large_certificate():
    # Independent reconstruction of one n=12 matching archived in round 82.
    n = 12
    seed = 821012
    edges = matching(n, seed)
    sha = hashlib.sha256(json.dumps(edges).encode()).hexdigest()
    N = 1 << (n - 1)
    parent_D = 3 * N + n + 1
    parent_s = 4 * N * N
    size = 1 << n
    masks = [None]
    for normal in range(1, size):
        b = [0, 0]
        for i, (x, y) in enumerate(edges):
            if (normal & (x ^ y)).bit_count() % 2 == 0:
                b[(normal & x).bit_count() % 2] |= 1 << i
        masks.append(tuple(b))

    def margin(a, f):
        M = size >> a
        t = M - 2 * f
        w = 2 - a
        D_upper = 3 * (f + t) + n - a + 1
        triple_s_lower = 6 * M * (f + 2 * t) - 16 * M * w
        return parent_D * parent_D * triple_s_lower - 3 * parent_s * D_upper * D_upper

    first = min(margin(1, b.bit_count()) for pair in masks[1:] for b in pair)
    second = None
    count = 0
    for a in range(1, size):
        for b in range(a + 1, size):
            if b > (a ^ b):
                continue
            for i in masks[a]:
                for j in masks[b]:
                    value = margin(2, (i & j).bit_count())
                    second = value if second is None or value < second else second
                    count += 1
    assert first > 0 and second > 0
    assert all(
        (4 * N * N - 2 * N * L) * parent_D**2 > parent_s * (parent_D - L)**2
        for L in range(1, 5)
    )

    labels = {x: i for i, edge in enumerate(edges) for x in edge}
    for a in range(3, n // 2 + 1):
        M = size >> a
        for answer in range(1 << a):
            cell = [(y << a) | answer for y in range(M)]
            active = {labels[y] for y in cell}
            f = M - len(active)
            t = M - 2 * f
            D = 3 * len(active) + n - a + 1
            s = 2 * M * (f + 2 * t)
            if t and D * D * parent_s > parent_D * parent_D * s:
                return sha, parent_D, parent_s, first, second, count, (a, answer, f, t, D, s)
    raise AssertionError("coordinate recovery not found")


def verify_local_and_sharp_cases():
    edges = [(0, 1), (2, 3), (4, 6), (5, 7)]
    cell = range(8)
    constants = direct_small_check(edges, cell, [(2, 0, 0)])
    split = direct_small_check(edges, cell, [(2, 1, 0)])
    identified = direct_small_check(edges, cell,
                                    [(9, 0, 0), (18, 0, 0), (36, 0, 0)])
    assert (constants, split, identified) == ((14, 52), (15, 64), (12, 40))
    return constants, split, identified


if __name__ == "__main__":
    print("Python", platform.python_version(), flush=True)
    print("internal constant, internal split, cross-bound equality:", verify_local_and_sharp_cases(), flush=True)
    print("n=3 direct one-query and sampled query counts:", verify_small_mixed(), flush=True)
    print("n=12 all h<=2, and explicit recovery:", large_certificate(), flush=True)
