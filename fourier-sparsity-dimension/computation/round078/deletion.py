"""Exhaust the four-bit Walsh supports and all linear cores."""
import json
from collections import Counter

def rank(columns):
    pivots = {}
    for x in columns:
        while x and x.bit_length() in pivots:
            x ^= pivots[x.bit_length()]
        if x:
            pivots[x.bit_length()] = x
    return len(pivots)

def span(gens):
    v = {0}
    for x in gens:
        v |= {t ^ x for t in v}
    return frozenset(v)

def spaces(n):
    known = {frozenset({0})}
    for _ in range(n):
        known |= {span([*h,x]) for h in known for x in range(1<<n) if x not in h}
    return known

def main():
    from functools import cache
    n=4
    chars=[sum(1<<x for x in range(16) if (x&a).bit_count()%2) for a in range(16)]
    cores=spaces(n)
    counts=Counter()
    @cache
    def eligible(support):
        ambient=span(support)
        tested=[]
        for h in cores:
            if not h<ambient: continue
            e=sorted(set(support)-h)
            if rank(16+x for x in e)==len(e):
                tested.append((tuple(e),rank(ambient)-rank(h)))
        return tested
    for table in range(1<<16):
        w=[16-2*(table^mask).bit_count() for mask in chars]
        s=tuple(i for i,z in enumerate(w) if z)
        for e,r in eligible(s):
            assert r==1 or (r==2 and len(e)==3 and len(s)==4)
            if len(s)>4:
                assert len(e)==2 and all(abs(w[i])==8 for i in e)
                u=e[0]^e[1]
                assert all((i^u) in s for i in s)
                # Fold the core modulo u and discard its redundant coordinate.
                base=u.bit_length()-1
                g={i^u if i>>base&1 else i for i in set(s)-set(e)}
                assert len(s)==2*len(g)+2 and rank(s)==rank(g)+2
            counts[f's={len(s)},r={r}']+=1
    print(json.dumps(dict(tables=65536,supports=eligible.cache_info().currsize,
                          subspaces=len(cores),checks=dict(sorted(counts.items()))),indent=2))
if __name__=='__main__':main()
