#!/usr/bin/env python3
"""Exact bounded arithmetic checks of R66 erosion and overlapping-charge claims."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import sys
import time


def main():
    start = time.perf_counter()
    scalar_count = 0
    for s in (F(1, 4), F(1, 2), F(1), F(2)):
        for k in range(5):
            t = k * s / 4
            for i in range(6):
                for j in range(6):
                    u, v = i * s / 4, j * s / 4
                    area = lambda q: max(q-u, 0) * max(q-v, 0)
                    diff = area(s) - area(t)
                    assert 0 <= diff <= s*s-t*t
                    scalar_count += 1
    families = []
    for n in (8, 16, 32, 64):
        xs = [F(3*i, 4*n) for i in range(n)] + [F(19, 16)]
        masses = [F(1, n)] * n + [F(1, 4)]
        eligible = []
        closure_count = 0
        for lo in range(n+1):
            for hi in range(lo, n+1):
                mass = sum(masses[lo:hi+1])
                d2 = (xs[hi]-xs[lo]) ** 2
                if d2 < mass:
                    eligible.append((lo, hi, mass))
                closure_count += 1
        L = max(row[2] for row in eligible)
        maxrows = [row[:2] for row in eligible if row[2] == L]
        assert L == 1 and maxrows == [(0, n-1)]
        assert (xs[-1]-xs[-2])**2 < masses[-1]+masses[-2]
        nonterminal = xs[-1]**2 - F(35, 33)*sum(masses)
        assert nonterminal > 0
        charges = []
        for k in range(n//4+1, n//2):
            a, b = F(n-k, n), F(1, 4)
            m, z = a+b, F(1, 4)-F(k, n)
            d = xs[-1]-xs[k]
            assert F(3, 4) < m < 1
            assert d == 1+3*z/4
            assert d*d-m == z/2+9*z*z/16 < 0
            assert a > F(1, 2)
            charges.append(a)
        total = sum(charges)
        assert total > F(n//4-1, 2)
        families.append({"N": n, "source_count": n+1,
                         "all_closed_intervals": closure_count,
                         "maximum_eligible_mass": "1", "unique_max_interval": [0,n-1],
                         "connected_by_full_S_and_last_pair": True,
                         "nonterminal_margin": str(nonterminal),
                         "crossing_high_closures": len(charges),
                         "sum_of_intersection_mass_capacities": str(total),
                         "strict_lower_for_capacity_sum": str(F(n//4-1, 2))})
    result = {"status": "PASS", "python": sys.version,
              "scalar_cases": scalar_count, "families": families,
              "arithmetic": "exact integer/Fraction; no randomness/floats",
              "scope": "Finite arithmetic supports analytic erosion bound and unbounded repeated-charge example; no general deletion theorem or halo-area calculation.",
              "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "seconds": time.perf_counter()-start}
    out=Path(sys.argv[1]);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({"status":"PASS","scalar_cases":scalar_count,"families":len(families),"seconds":result["seconds"]}))


if __name__ == '__main__':
    main()
