"""Check the first-merge reflection inequalities using exact suffix counts."""
from collections import Counter
from functools import cache
from math import factorial
import json
import platform
import random


def check(n, howmany):
    N = n*n
    full = (1 << N)-1
    mid = n//2
    w = mid*n+mid
    seam = sum(1 << (y*n+mid) for y in range(n))
    left = sum(1 << v for v in range(N) if v%n < mid)
    right = full ^ left ^ seam
    neighbors = [0]*N
    for v in range(N):
        for u in range(N):
            if abs(u%n-v%n)+abs(u//n-v//n)==1:
                neighbors[v] |= 1 << u

    @cache
    def around(mask):
        if not mask:
            return 0
        low = mask & -mask
        return around(mask^low) | neighbors[low.bit_length()-1]

    @cache
    def suffix(mask):
        if mask == full:
            return 1
        return sum(suffix(mask|b) for b in bits(around(mask)&~mask))

    def reflected(mask):
        return sum(1 << ((n-1-v//n)*n+v%n) for v in vertices(mask))

    def random_component(side, root, remove):
        out = side
        while remove:
            u = rng.choice(tuple(v for v in vertices(out) if v != root))
            smaller = out^(1<<u)
            component = 1<<root
            while True:
                nxt = component | (around(component)&smaller)
                if component == nxt:
                    break
                component = nxt
            if component==smaller:
                out = smaller
                remove -= 1
        return out

    rng = random.Random(49182026+n)
    cases = Counter()
    for case in range(howmany):
        a = random_component(left,w-1,rng.randrange(max(1,n-1)))
        b = random_component(right,w+1,rng.randrange(max(1,n-1)))
        if case < 2:
            a,b = left,right
        s, t = a|b|(1<<w), a|reflected(b)|(1<<w)
        assert t.bit_count()==s.bit_count()
        cleared = suffix(s|seam)
        assert cleared == suffix(t|seam)
        m, k = N-s.bit_count(), n-1
        ways = factorial(m)//factorial(m-k)
        x,y = suffix(s),suffix(t)
        assert cleared <= x <= cleared*ways
        assert cleared <= y <= cleared*ways
        assert x <= ways*y and y <= ways*x
        for first, c in ((s,x),(t,y)):
            if (seam&~first)&~around(first)==0:
                assert c >= factorial(k)*cleared
        u = full&~(s|around(s))
        v = full&~(t|around(t))
        band = seam | around(seam)
        d = v&band
        h = v&~band
        if not (u|v)&band:
            assert x==y
            cases['exact'] += 1
        if not around(d)&(d|h):
            assert y*3**d.bit_count() >= x*2**d.bit_count()
            cases['sparse'] += 1
        cases['checked'] += 1
        if x!=y:
            cases['unequal'] += 1
    return dict(n=n,cases=cases,suffix_states=suffix.cache_info().currsize)


def bits(mask):
    while mask:
        bit=mask&-mask
        yield bit
        mask^=bit


def vertices(mask):
    return [x.bit_length()-1 for x in bits(mask)]


if __name__=='__main__':
    print(json.dumps(dict(python=platform.python_version(), checks=[check(3,100),check(5,120),check(7,60)]),indent=2))
