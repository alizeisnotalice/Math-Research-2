from collections import Counter
from fractions import Fraction
from itertools import product
from random import Random
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


def support(values):
    return tuple(i for i, value in enumerate(values) if value)


def pivots(vectors):
    result = {}
    for x in vectors:
        while x:
            key = x.bit_length()
            if key not in result:
                result[key] = x
                break
            x ^= result[key]
    return result


def rank(vectors):
    return len(pivots(vectors))


def reduce_mod(x, basis):
    for key, v in sorted(basis.items(), reverse=True):
        if x >> (key - 1) & 1:
            x ^= v
    return x


def affine_rank(points, dim):
    return rank((1 << dim) | point for point in points)


def mask_gap(values):
    N = len(values)
    active = support(values)
    if not active:
        return None
    physical_size = len(active)
    spectrum = support(walsh(values))
    product_count = physical_size * len(spectrum)
    if product_count == N:
        d = N.bit_length() - 1
        assert physical_size == 1 << (affine_rank(active, d) - 1)
        phase = 1 if values[active[0]] > 0 else -1
        assert any(all(values[x] == phase * (1 - 2*((w & (x ^ active[0])).bit_count() % 2))
                       for x in active) for w in range(N))
    else:
        assert 4 * product_count >= 9 * N
    return Fraction(product_count, N)


def verify_gap():
    exact = Counter(mask_gap(values) for values in product((-1,0,1), repeat=8))
    rng = Random(8601)
    sampled = 0
    for _ in range(5000):
        d = 4 + rng.randrange(2)
        values = [rng.randrange(3)-1 for _ in range(1 << d)]
        mask_gap(values)
        sampled += 1
    sharp = [int(x & 3 == 0) - int(x >> 2 == 0) for x in range(16)]
    assert mask_gap(sharp) == Fraction(9, 4)
    return sum(exact.values()) - exact[None], sampled, str(mask_gap(sharp))


def check_node(naddr, roles, signs, children):
    if len(roles) != 1 << naddr:
        raise AssertionError("addresses not fully partitioned")
    labels = sorted(set(roles))
    Ns = [len(child) for child in children]
    ds = [s.bit_length()-1 for s in Ns]
    data_support = [support(walsh(values)) for values in children]
    ss = [len(a) for a in data_support]
    Ri = [affine_rank(a, dim) for a, dim in zip(data_support, ds)]
    offsets = []
    d = naddr
    for bits in ds:
        offsets.append(d)
        d += 1+bits
    masks = {e: [sign if branch == e else 0 for branch, sign in zip(roles, signs)]
             for e in labels}
    bands = {e: support(walsh(mask)) for e, mask in masks.items()}
    W = pivots(anchor ^ band[0] for band in bands.values() for anchor in band)
    lifted = rank((1 << (naddr+len(children))) | (e << naddr) | reduce_mod(band[0], W)
                  for e, band in bands.items())
    expected_R = len(W) + lifted + sum(R-1 for R in Ri)
    expected_s = 0
    for e, band in bands.items():
        width = len(band)
        for i, sz in enumerate(ss):
            if (e >> i) & 1:
                width *= sz
        expected_s += width
    values = []
    for point in range(1 << d):
        addr = point & ((1 << naddr)-1)
        e = roles[addr]
        sign = signs[addr]
        for i, child in enumerate(children):
            if (e >> i) & 1:
                x = point >> (1+offsets[i]) & ((1 << ds[i])-1)
                guard = point >> offsets[i] & 1
                sign *= child[x] * (1-2*guard)
        values.append(sign)
    S = support(walsh(values))
    R = affine_rank(S, d)
    assert (len(S), R) == (expected_s, expected_R)
    m = len(labels)
    base = 1 << (m.bit_length()-1)
    delta = sum(map(len, bands.values()))-base*(3*m-2*base)
    assert len(W) <= m-1+2*delta
    assert sum(len(bands[e]) for e in labels if mask_gap(masks[e]) > 1) <= 2*delta
    T = Fraction(1 << naddr, m)
    tau = 4*T-3
    K2 = max((4*T-2)**2, *[Fraction((r+tau)**2, s) for r,s in zip(Ri,ss)])
    assert Fraction((R+tau)**2, len(S)) <= K2
    return values, R, len(S), delta


def verify_graphs():
    rng = Random(8602)
    seeds = ((1,), (1,-1), (1,1,1,-1))
    configs = [
        (2, (3,3,3,3), (1,1,1,-1), ((1,-1),(1,1,1,-1))),
        (2, (0,1,3,3), (1,-1,-1,1), ((1,1,1,-1),(1,-1))),
    ]
    for _ in range(28):
        children = tuple(rng.choice(seeds) for _ in range(3))
        edges = [3, 5, 6, 7]
        roles = [rng.choice(edges) for _ in range(8)]
        for index, edge in enumerate(roles):
            if not edge:
                roles[index] = 7
        signs = [rng.choice((-1,1)) for _ in roles]
        configs.append((3,roles,signs,children))
    outcomes = [check_node(*config) for config in configs]
    grandchild = check_node(2, (1,1,3,3), (1,1,1,-1), (outcomes[0][0], (1,)))
    return len(outcomes), (grandchild[1], grandchild[2], grandchild[3])


def two_edge_encoder():
    n = 6
    rng = Random(8603)
    values = [rng.choice((-1,1)) for _ in range(1 << n)]
    W = walsh(values)
    d = rank(support(W))
    roles = [3 if val == 1 else 5 for val in values]
    _, R, s, _ = check_node(n, roles, [1]*len(roles), [(1,), (1,), (1,)])
    assert R == d+2
    assert s == 2*(len(support(W))+(W[0] == 0))
    return d, len(support(W)), R, s


if __name__ == "__main__":
    print("Python", platform.python_version(), flush=True)
    print("mask gap: exhaustive nonzero ternary tables, random profiles, sharp value:", verify_gap(), flush=True)
    print("actual guarded graph/hypergraph nodes and one two-level result:", verify_graphs(), flush=True)
    print("unsigned two-edge encoder (d_G,s_G,R_H,s_H):", two_edge_encoder(), flush=True)
