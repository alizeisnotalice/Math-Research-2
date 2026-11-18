"""Run python3 -B check_weak_folding.py. All identities use integers."""
import itertools
import json
import platform
from collections import Counter, defaultdict
from fractions import Fraction

MODULI = {2: 0b111, 3: 0b1011, 4: 0b10011, 5: 0b100101}


def mul(x, y, q, modulus):
    result = 0
    while y:
        if y & 1:
            result ^= x
        x <<= 1
        y >>= 1
        if x & q:
            x ^= modulus
    return result


def independent(columns):
    pivots = {}
    retained = []
    for x in columns:
        v = x
        while v and v.bit_length() in pivots:
            v ^= pivots[v.bit_length()]
        if v:
            pivots[v.bit_length()] = v
            retained.append(x)
    return retained


def make_blocks(k):
    q = 2**k
    modulus = MODULI[k]
    def trace(x):
        r = 0
        for _ in range(k):
            r ^= x
            x = mul(x, x, q, modulus)
        assert r in (0, 1)
        return r
    a = next(x for x in range(1, q) if trace(x))
    x0 = next(x for x in range(q) if trace(mul(a, x, q, modulus)))
    blocks = []
    for i in range(q):
        tag = 1 << (2*k+i)
        block = {}
        for b in range(q):
            c = mul(i, b, q, modulus)
            address = min(c, c ^ a) | (b << k)
            block[tag | address] = (-1)**trace(mul(c, x0, q, modulus))
        blocks.append(block)
    return blocks


def grow(subspace, v):
    return frozenset(subspace | {x ^ v for x in subspace})


def audit(k, exhaustive):
    q = 2**k
    blocks = make_blocks(k)
    full = {x: sign for block in blocks for x, sign in block.items()}
    rank = len(independent(full))
    assert len(full) == q*q and rank == q+2*k-1
    signed_sums = Counter()
    for (x, a), (y, b) in itertools.combinations(full.items(), 2):
        signed_sums[x ^ y] += a*b
    assert not any(signed_sums.values())
    # Physical address masks give all characters on each queried kernel.
    mask = (1 << (2*k))-1
    address_space = {0}
    for v in independent(x & mask for x in full):
        address_space = grow(address_space, v)
    assert len(address_space) == q*q//2
    queue = {frozenset({0})}
    snapshots = []
    terminal = None
    for depth in range(2*k):
        next_queue = set()
        local_stats = []
        for kernel in queue:
            projection = {a: min(a ^ h for h in kernel) for a in address_space}
            fibers = defaultdict(list)
            for x, sign in full.items():
                v = (x & ~mask) | projection[x & mask]
                fibers[v].append((x & mask, sign))
            projected = sorted(fibers)
            fold = Counter(x ^ y for x, y in itertools.combinations(projected, 2))
            n = sorted((len({projection[x & mask] for x in block}) for block in blocks),
                       reverse=True)
            w = 2*k-1-depth
            r1 = n[0].bit_length()-1
            expected = Fraction(sum(v*(v-1) for v in n)-n[0]*r1-n[1]*(w-r1), q)
            # Independently transform each signed fiber over kernel assignments.
            basis_h = independent(sorted(kernel))
            characters = {}
            for physical in range(1 << (2*k)):
                key = tuple((physical & v).bit_count() % 2 for v in basis_h)
                characters.setdefault(key, physical)
            hits = {}
            for v, fiber in fibers.items():
                hits[v] = sum(bool(sum(sign*(-1)**((a & x).bit_count() % 2)
                                       for a, sign in fiber)) for x in characters.values())
            basis = independent(sorted(hits, key=hits.get, reverse=True))
            observed = Fraction(sum(hits.values())-sum(hits[x] for x in basis), len(kernel))
            assert observed == expected
            if not w:
                assert expected == 0
                terminal = (depth, len(fibers))
                continue
            address_best = max(value for direction, value in fold.items() if direction <= mask)
            mixed_best = max((value for direction, value in fold.items() if direction > mask),
                             default=0)
            assert address_best > mixed_best
            local_stats.append((address_best-mixed_best, expected))
            if exhaustive:
                candidates = set(projection.values())-{0}
            else:
                candidates = [min(v for v, value in fold.items()
                                  if v <= mask and value == address_best)]
            for direction in candidates:
                next_queue.add(grow(kernel, direction))
        snapshots.append({"depth": depth, "subspaces": len(queue),
                          "min_greedy_margin": min((margin for margin, _ in local_stats), default=None),
                          "min_Et": str(min((et for _, et in local_stats), default=0)),
                          "max_Et": str(max((et for _, et in local_stats), default=0))})
        queue = next_queue
    return {"k": k, "Q": q, "s": q*q, "d": rank,
            "c": (q*q+rank-1)//rank, "scope": "all address subspaces" if exhaustive else "one greedy path",
            "signed_nonzero_differences_checked": len(signed_sums), "levels": snapshots,
            "terminal_depth_and_cosets": terminal}


if __name__ == "__main__":
    print(json.dumps({"python": platform.python_version(), "seed": None,
                      "tolerance": 0, "instances": [audit(k, k <= 3) for k in (2, 3, 4, 5)]},
                     indent=2))
