"""Exact contractions and greedy bases for the short-sumset bound."""
from collections import defaultdict
from fractions import Fraction
from math import comb, isqrt
import json
import platform


def minimum_count(j, k):
    def ball(m):
        return sum(comb(m, i) for i in range(min(k, m) + 1))
    lo, hi = 0, 1
    while ball(hi) < 1 << j:
        hi *= 2
    while hi-lo > 1:
        mid = (lo+hi)//2
        if ball(mid) >= 1 << j:
            hi = mid
        else:
            lo = mid
    return hi


def shifted(mask, v):
    return sum(1 << (x^v) for x in range(16) if mask & (1 << x))


def enlarge(mask, v):
    return mask | shifted(mask, v)


def lattice():
    spans = {1}
    for x in range(1, 16):
        spans |= {enlarge(h,x) for h in spans}
    return spans


def sums(labels, k):
    layers = [1] * (k+1)
    for b in labels:
        for j in range(k, 0, -1):
            layers[j] |= shifted(layers[j-1], b)
    return layers[k]


def exact_greedy(labels):
    states = {1: Fraction(1)}
    marginals = defaultdict(Fraction)
    while states:
        following = defaultdict(Fraction)
        for h, mass in states.items():
            available = [e for e in labels if not h >> e & 1]
            for e in available:
                prob = mass / len(available)
                marginals[e] += prob
                following[enlarge(h,e)] += prob
        states = following
    return sum(x*x for x in marginals.values())


def main():
    subs = lattice()
    checks = 0
    # Every nonzero label subset in four dimensions. Pair with every
    # contraction whose kernel lies in its span. Use k=2 to test rank loss.
    floors = [0] + [minimum_count(j,2) for j in range(1,5)]
    for mask in range(0, 1 << 16, 2):
        labels = [e for e in range(1,16) if mask >> e & 1]
        h = 1
        for e in labels:
            h = enlarge(h,e)
        d = h.bit_count().bit_length()-1
        reached = sums(labels,2).bit_count()
        p = d - (reached.bit_length()-1)
        for sub in subs:
            if sub & ~h:
                continue
            j = d - (sub.bit_count().bit_length()-1)
            if j > p:
                assert (mask & ~sub).bit_count() >= floors[j-p]
            checks += 1
    scale = 10**8
    upper = 0
    for j in range(1,129):
        m = minimum_count(j,4)
        upper += isqrt((scale*scale)//m) + 1
    # m_4(j) >= 2^(j/4)-1, with 2^(-1/8)<23/25.
    assert 2 * 23**8 > 25**8
    upper4 = Fraction(upper,scale) + Fraction(25,65536)
    samples = [[1,2,4,8], list(range(1,16)), [1,2,3,4,8], [1,2,4,8,15]]
    return dict(python=platform.python_version(), labels_tested=32768,
                contractions_tested=checks, exact_upper_A4=str(upper4),
                upper_A4_rounded_up=(int(upper4*1000)+1)/1000,
                greedy_squared_costs=[str(exact_greedy(s)) for s in samples])


if __name__ == '__main__':
    print(json.dumps(main(), indent=2))
