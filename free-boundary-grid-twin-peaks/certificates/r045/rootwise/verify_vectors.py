"""Reconstruct exact root-stratified lower bounds from sparse root vectors."""
import json
from collections import Counter, defaultdict
from fractions import Fraction
from math import comb, factorial
from pathlib import Path

here=Path(__file__).resolve().parent
for n in (7,8):
    vectors=json.loads((here/f'n{n}-root-vectors.json').read_text())
    stored=json.loads((here/f'n{n}-summary.json').read_text())
    def root_vector(h):
        if h in vectors:
            return [(p,int(a)) for p,a in vectors[h]]
        return [((n-1-p//n)*n+p%n,int(a)) for p,a in vectors[h[::-1]]]
    def pair_after(pair,f):
        return tuple(sorted(f(v%n,v//n) for v in pair))
    if n==7:
        family=stored['event_counts']
    else:
        family={z['profile']:z['count'] for z in stored['event_counts']}
    counts=defaultdict(Counter)
    for h,expected in family.items():
        parity=tuple(int(w)%2 for w in h)
        reflected=''.join(str(n-int(w)) for w in h)
        shuffle=comb(n*n,sum(map(int,h)))
        total=0
        for left,a in root_vector(h):
            for right,b in root_vector(reflected):
                rx,ry=n-1-right%n,right//n
                assert abs(rx-left%n)+abs(ry-left//n)>n/2
                pair=tuple(sorted((left,n*ry+rx)))
                value=shuffle*a*b
                counts[parity,0][pair]+=value
                total+=value
        assert total==int(expected)
    originals=list(counts.items())
    for (tag,_),mass in originals:
        counts[tag,1]=Counter({pair_after(r,lambda x,y:n*x+y): v for r,v in mass.items()})
        if n==8:
            reverse=tuple(reversed(tag))
            counts[reverse,0]=Counter({pair_after(r,lambda x,y:n*(n-1-y)+x):v for r,v in mass.items()})
            counts[reverse,1]=Counter({pair_after(r,lambda x,y:n*x+n-1-y):v for r,v in mass.items()})
    maxima={}
    for mass in counts.values():
        for pair,v in mass.items():
            maxima[pair]=max(v,maxima.get(pair,0))
    total=sum(maxima.values())
    if n==7:
        assert total==int(stored['rootwise_disjoint_count'])
        assert {tuple(d['peaks']):int(d['count']) for d in stored['root_strata']}==maxima
        assert Fraction(total,int(stored['A2']))==Fraction(stored['rootwise_disjoint_probability_exact'])
    else:
        assert total==int(stored['rootwise_disjoint_subunion_count'])
        assert {tuple(d['peaks']):int(d['chosen_count']) for d in stored['root_strata']}==maxima
        assert Fraction(total,factorial(n*n))==Fraction(stored['unconditional_rootwise_probability_exact'])
    print(n,len(vectors),len(family),len(maxima),total)
