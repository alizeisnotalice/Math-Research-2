"""python3 -B check_unguarded.py. Exhaust all two-bit two-slot signed masks."""
import json
import platform
from collections import Counter
from itertools import product
from fractions import Fraction


def walsh(values):
    out = list(values)
    step = 1
    while step < len(out):
        for start in range(0,len(out),2*step):
            for j in range(start,start+step):
                a,b = out[j],out[j+step]
                out[j],out[j+step] = a+b,a-b
        step *= 2
    return out


def rank(columns):
    pivots = {}
    for x in columns:
        while x and x.bit_length() in pivots:
            x ^= pivots[x.bit_length()]
        if x:
            pivots[x.bit_length()] = x
    return len(pivots)


def stats(values):
    S = [i for i,v in enumerate(walsh(values)) if v]
    return len(S)-(0 in S),1+rank(2*x+1 for x in S)


def tables():
    return [[1-2*((bits>>i)&1) for i in range(4)] for bits in range(16)]


def main():
    a = 2
    c = 3*(a+2)
    rows = tables()
    childstats = [stats(row) for row in rows]
    checks = Counter()
    for slots in product((0,1),repeat=4):
        # A single child left unconstrained below accounts for unary masks too.
        labels = set(slots)
        if labels == {1}:
            continue
        for signs in rows:
            bands = [walsh(signs[z] if slots[z]==i else 0 for z in range(4))
                     for i in labels]
            bs = [sum(x!=0 for x in band) for band in bands]
            assert sum(Fraction(1,b) for b in bs)<=1
            for seeds in product(range(16),repeat=len(labels)):
                stats_i = [childstats[s] for s in seeds]
                active = [(k,D) for k,D in stats_i if k]
                # Compute Walsh spectrum directly from the physical truth table.
                vals = [signs[z]*rows[seeds[slots[z]]][(inputs>>(2*slots[z]))&3]
                        for inputs in range(4**len(seeds)) for z in range(4)]
                k,D = stats(vals)
                checks[f'active_{len(active)}']+=1
                assert k>=sum(b*ki for b,(ki,_) in zip(bs,stats_i))
                if not k:
                    assert not active
                    continue
                if active:
                    assert any((D+c)**2*ki<=(Di+c)**2*k for ki,Di in active)
                else:
                    assert D<=a+2
                if all(sum(rows[s]) for s in seeds):
                    mean_address = [signs[z]*sum(rows[seeds[slots[z]]]) for z in range(4)]
                    if 0 in walsh(mean_address):
                        checks['all_biased_some_zero_band_cancellation']+=1
    print(json.dumps({'python':platform.python_version(),'seed':None,'tolerance':0,
                      'scope':'all signed assignments of four addresses to one or two labeled nonempty slots, all two-bit seeds',
                      'checks':dict(checks)},indent=2))


if __name__=='__main__':
    main()
