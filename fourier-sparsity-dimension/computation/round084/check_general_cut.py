"""True Walsh bands for quadratic pencils and all cuts on small parents."""
import json
import platform
import random
from fractions import Fraction
from cut_mass import cost


def powers(bits):
    return [1 << j for j in range(bits)]


def one(seed_dim, address_dim, pencils):
    monomials = [a | b for a in powers(address_dim) for b in powers(address_dim) if a < b]
    seed_values = ([1] if not seed_dim else [1, 1, 1, -1])
    data_count, addr_count = 1 << seed_dim, 1 << address_dim
    child, parent, ranks = {}, {}, set()
    for e in range(data_count):
        old = abs(sum(value * (-1)**((x & e).bit_count() % 2) for x, value in enumerate(seed_values)))
        v = data_count | e
        child[v] = old
        polynomial = 0
        for i, p in enumerate(pencils):
            if v & (1 << i):
                polynomial ^= p
        val = [(-1)**(sum(1 for j, m in enumerate(monomials)
                            if polynomial >> j & 1 and x & m == m) % 2) for x in range(addr_count)]
        band = [abs(sum(value * (-1)**((x & a).bit_count() % 2) for x, value in enumerate(val)))
                for a in range(addr_count)]
        ranks.add(sum(bool(w) for w in band).bit_length() - 1)
        for a, w in enumerate(band):
            if w:
                parent[v | (a << (seed_dim + 1))] = w * old
    old_cost, c = cost(child)
    fresh, q = cost(parent)
    ratio = addr_count * fresh / old_cost
    bound = Fraction(3, 2) if 2 in ranks else Fraction(5, 4) if 4 in ranks else 1
    assert ratio <= bound
    if address_dim <= 3 or ranks <= {0, 2}:
        # Full affine label domain consists of exactly the tested seed spectrum.
        assert ratio <= 1
    return q + c


if __name__ == '__main__':
    rng = random.Random(840801)
    records = []
    for n, m, trials in ((2, 3, 96), (2, 4, 96), (0, 6, 24)):
        checks = sum(one(n, m, [rng.getrandbits(m * (m - 1) // 2) for _ in range(n + 1)])
                     for _ in range(trials))
        records.append(dict(seed_bits=n, address_bits=m, cases=trials, checked_subspaces=checks))
    print(json.dumps(dict(python=platform.python_version(), seed=840801, records=records), indent=2))
