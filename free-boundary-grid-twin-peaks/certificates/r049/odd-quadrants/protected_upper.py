"""Four-block necessary-condition denominator and the odd-grid lower bound."""
from fractions import Fraction
from functools import cache
from itertools import permutations, product
from math import comb, factorial, prod
from pathlib import Path
import json
import platform
from verify import bits

P = Path(__file__).resolve().parent

def side(height, width):
    n = height*width
    full = (1<<n)-1
    adj = [sum(1<<w for w in range(n)
               if abs(v%width-w%width)+abs(v//width-w//width)==1) for v in range(n)]
    protected = sum(1 << (y*width+x) for y in range(height-1) for x in range(width-1))
    outside = full ^ protected

    @cache
    def closed(S):
        D = S
        for b in bits(S):
            D |= adj[b.bit_length()-1]
        return D

    @cache
    def simultaneous(S):
        if not S:
            return factorial(n)
        value, remainder = divmod(sum(simultaneous(S^b) for b in bits(S)),
                                  closed(S).bit_count())
        assert not remainder
        return value

    moments = [[0]*(protected.bit_count()+1) for _ in range(3)]

    def visit(S, left):
        k = S.bit_count()
        D = closed(S)
        a, b = protected&~D, outside&~D
        ea = sum((adj[j.bit_length()-1]&b).bit_count() for j in bits(a))
        base = simultaneous(S)
        moments[0][k] += base
        t,rt = divmod(b.bit_count()*base,n)
        typ = ((a.bit_count()+D.bit_count())*b.bit_count()-ea)*base
        v,rv = divmod(typ,n*(n-1))
        assert not rt and not rv
        moments[1][k] += t
        moments[2][k] += v
        while left:
            chosen=left&-left
            left ^= chosen
            visit(S|chosen,left&~adj[chosen.bit_length()-1])

    visit(0,protected)
    histograms = [[sum((-1)**(k-j)*comb(k,j)*m[k] for k in range(j,len(m)))
                   for j in range(len(m))] for m in moments]
    f,t,v=histograms
    assert sum(f)==factorial(n)
    assert sum(t)==outside.bit_count()*factorial(n-1)
    assert all(0 <= x <= u-y <= u for x,u,y in zip(v,f,t,strict=True))
    return histograms, simultaneous.cache_info().currsize

def brute(height,width):
    n = width*height
    adj = [[w for w in range(n)
            if abs(w%width-v%width)+abs(w//width-v//width)==1] for v in range(n)]
    eligible = {y*width+x for y in range(height-1) for x in range(width-1)}
    f,t,v=([0]*(len(eligible)+1) for _ in range(3))
    for w in permutations(range(n)):
        rank={item:i for i,item in enumerate(w)}
        i=sum(all(rank[e]<rank[u] for u in adj[e]) for e in eligible)
        f[i]+=1
        t[i]+=w[0] not in eligible
        v[i]+=w[0] in eligible and w[1] not in eligible and w[1] not in adj[w[0]]
    return [f,t,v]

def multinomial(areas):
    return factorial(sum(areas))//prod(factorial(a) for a in areas)

def count_upper(areas, components):
    names = range(len(areas))
    def prefix(*fixed):
        a=list(areas)
        for i in fixed:
            a[i]-=1
        return multinomial(a)

    raw = 0
    first = 0
    second = 0
    for ks in product(range(3),repeat=len(areas)):
        if sum(ks)>2:
            continue
        f,t,v=([components[i][j][ks[i]] for i in names] for j in range(3))
        raw += prefix()*prod(f)
        if sum(ks)==2:
            for j in names:
                other=prod(f[i] for i in names if i!=j)
                first += prefix(j)*t[j]*other
                second += prefix(j,j)*v[j]*other
                for k in names:
                    if k!=j:
                        second += prefix(j,k)*(f[j]-t[j])*t[k]*prod(f[i] for i in names if i not in (j,k))
    assert 0 < raw-first-second <= raw <= factorial(sum(areas))
    return raw, first, second, raw-first-second

def main():
    for x in ((2,2),(2,3),(3,3)):
        assert side(*x)[0]==brute(*x)
    shapes=[(5,5),(6,5),(5,6),(6,6)]
    sides=[side(*s) for s in shapes]
    assert sides[1]==sides[2]
    areas=[h*w for h,w in shapes]
    before, exclude_first, exclude_second, upper=count_upper(areas,[x[0] for x in sides])
    numerator=json.loads((P/'n11-numerator.json').read_text())['constructed_words']
    conditional=Fraction(numerator,upper)
    answer=dict(python=platform.python_version(),n=11,both_protected_peaks_upper=upper,
                initial=before,exclude_global_max=exclude_first,exclude_global_second=exclude_second,
                upper_probability=str(Fraction(upper,factorial(sum(areas)))),
                known_lower_far_threshold='11/2',conditional_lower=str(conditional),
                decimal=float(conditional),brute_side_sizes=['2x2','2x3','3x3'],
                regional_independent_sets=[x[1] for x in sides],sides=[x[0] for x in sides])
    print(json.dumps(answer,indent=2))

if __name__=='__main__':
    main()
