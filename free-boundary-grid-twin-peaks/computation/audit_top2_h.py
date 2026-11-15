"""Independent Python-int audit of h({u,v}) for n<=5.

It evaluates h by the global-top-two inclusion-exclusion identity from the
joint-peak counts C(I), rather than by the reverse growth DP in top2_h_dp.
"""
from collections import defaultdict
from functools import cache
from math import factorial
import json
from pathlib import Path
import subprocess
import sys


def independent_sets(n):
    rows = [r for r in range(1 << n) if not (r & (r << 1))]
    answer = []

    def visit(y, old, whole):
        if y == n:
            answer.append(whole)
            return
        for row in rows:
            if not row & old:
                visit(y + 1, row, whole | (row << (y*n)))
    visit(0, 0, 0)
    return answer


def audit(n, binary):
    N = n*n
    full = (1 << N) - 1

    def nb(s):
        z = 0
        for v in range(N):
            if s >> v & 1:
                x, y = v % n, v // n
                for X, Y in ((x-1,y), (x+1,y), (x,y-1), (x,y+1)):
                    if 0 <= X < n and 0 <= Y < n:
                        z |= 1 << (Y*n+X)
        return z

    sets = independent_sets(n)

    @cache
    def c(s):
        if not s:
            return factorial(N)
        numerator = sum(c(s ^ (1 << v)) for v in range(N) if s >> v & 1)
        den = s.bit_count() + nb(s).bit_count()
        q, r = divmod(numerator, den)
        assert not r
        return q

    values = {s: c(s) for s in sets}
    bydist = defaultdict(int)
    for u in range(N):
        for v in range(u + 1, N):
            d = abs(u % n-v % n) + abs(u // n-v // n)
            if d == 1:
                continue
            seed = (1 << u) | (1 << v)
            z = full ^ (seed | nb(seed))
            numerator = sum((-1 if s.bit_count() & 1 else 1)*q
                            for s, q in values.items() if not s & ~z)
            q, r = divmod(numerator, N*(N-1))
            assert not r and q > 0
            bydist[d] += q
    result = {'n': n, 'ordered_top2_sum': str(2*sum(bydist.values())),
              'by_distance': {str(d): {'h_sum': str(q)} for d, q in bydist.items()}}
    cpp = json.loads(subprocess.check_output([str(binary), str(n)], text=True))
    assert result['ordered_top2_sum'] == cpp['ordered_top2_sum']
    assert result['by_distance'] == {d: {'h_sum': row['h_sum']}
                                     for d, row in cpp['by_distance'].items()}
    return result


def main(binary):
    out = [audit(n, binary) for n in (3, 4, 5)]
    print(json.dumps({'status': 'passed', 'results': out}, indent=2))


if __name__ == '__main__':
    main(Path(sys.argv[1]))
