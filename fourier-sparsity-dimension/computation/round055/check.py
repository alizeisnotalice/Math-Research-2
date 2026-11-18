"""Exact complete truth-table checks for ordinary full-coordinate spaces."""
from collections import Counter
from fractions import Fraction
import json
import platform

def rank(xs):
    pivots = {}
    for x in xs:
        while x and x.bit_length() in pivots:
            x ^= pivots[x.bit_length()]
        if x:
            pivots[x.bit_length()] = x
    return len(pivots)

def census(n):
    N = 1 << n
    masks = [sum(1 << x for x in range(N) if (a & x).bit_count() % 2) for a in range(N)]
    windows = [[a for a in range(1,N) if (v & a).bit_count() % 2] for v in range(1,N)]
    histogram = Counter()
    faithful = 0
    for f in range(1 << N):
        w = [N-2*(f ^ m).bit_count() for m in masks]
        if rank(a for a in range(N) if w[a]) != n:
            continue
        faithful += 1
        total = sum(abs(i) for i in w[1:])
        deficit = total - (2*N-2)
        assert deficit >= 0
        histogram[deficit] += 1
        if deficit == 0:
            assert all(abs(i) == 2 for i in w[1:])
            assert f.bit_count() in {1,N-1}
        else:
            assert deficit >= N//4 if n == 3 else deficit >= N//2
        if n == 4 and deficit:
            bad,good = [],[]
            EE = []
            for v,win in enumerate(windows,1):
                e = sum(abs(w[a]) for a in win)-N
                EE.append(e)
                (bad if e else good).append(v)
                assert e == 0 or 2*e >= N
            assert rank(bad)==n
            if rank(good)==n:
                assert all(e % N == 0 for e in EE)
                if sum(EE)>N:
                    assert 3*sum(EE)>=N*(N-1)
    return {"rank":n,"truth_tables":1<<N,"faithful":faithful,"endpoint":histogram[0],
            "least_positive_gap":str(Fraction(min((k for k in histogram if k),default=0),N)) if len(histogram)>1 else None,
            "distribution":{str(Fraction(k,N)):v for k,v in sorted(histogram.items())}}

print(json.dumps({"python":platform.python_version(),"normalization":"Walsh/N","ordinary_full_physical_space_only":True,
                  "census":[census(n) for n in range(1,5)]},indent=2))
