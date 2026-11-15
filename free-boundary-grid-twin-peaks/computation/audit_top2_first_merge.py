"""Check the first-merge interleaving formula for the global-top-two sector."""
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent))
from exact import component, grid, tables, vertices


EXPECTED = {3: 81984, 4: 350682265920}


def count(n):
    xy, adj, one, two, h, roots = tables(n, n, root_vectors=False)
    total = 0
    for S in range(1, 1 << (n*n)):
        if not one[S] or S.bit_count() < 3:
            continue
        for w in vertices(S):
            T = S ^ (1 << w)
            A = component(T, adj)
            B = T ^ A
            if not B or component(B, adj) != B:
                continue
            a, b = A.bit_count(), B.bit_count()
            # First two interleaving positions must be the roots of A,B.
            total += h[S]*one[A]*one[B]*2*math.comb(a+b-2, a-1)
    return total


for n, expected in EXPECTED.items():
    got = count(n)
    assert got == expected
    print({'n': n, 'top2_from_first_merge': got})
