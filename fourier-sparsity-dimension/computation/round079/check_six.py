"""Exact ternary-circuit and Boolean completion verification, no dependencies."""
import json
import platform
from collections import Counter, defaultdict
from itertools import combinations, product


def walsh(values):
    v = list(values)
    k = 1
    while k < len(v):
        for a in range(0, len(v), 2*k):
            for i in range(a, a+k):
                u, t = v[i], v[i+k]
                v[i], v[i+k] = u+t, u-t
        k *= 2
    return v


def rank(values):
    pivots = {}
    for t in values:
        while t and t.bit_length() in pivots:
            t ^= pivots[t.bit_length()]
        if t:
            pivots[t.bit_length()] = t
    return len(pivots)


def small():
    fixed = {}
    for n, target in [(4, {0,1,2,3}), (8, {0,1,2,3,4})]:
        counts = Counter()
        for f in product((-1,0,1), repeat=n):
            w = walsh(f)
            if {a for a, b in enumerate(w) if b} != target:
                continue
            z = [i for i in range(n) if not f[i]]
            counts[len(f)-len(z)] += 1
            if n == 8:
                for signs in product((-1,1),repeat=len(z)):
                    h = [0]*n
                    for i, t in zip(z, signs):
                        h[i] = t
                    h_spec = walsh(h)
                    assert sorted(abs(b) for b in h_spec if b) == [2,2,2,2,4]
                    for extra in (0,1):
                        S = [a+n for a in target]+[
                            a+extra*2*n for a, b in enumerate(h_spec) if b]
                        assert len(S) == 10 and rank(S) == 4+extra
        fixed[n] = counts
    return fixed


def six():
    X = (0,1,2,4,8,15)
    totals = Counter()
    model = None
    for coeff in product((-1,1),repeat=6):
        if sum(c<0 for c in coeff)%2 == 0:
            continue
        base = [0]*16
        for a, c in zip(X,coeff):
            base[a] = c
        raw = walsh(base)
        assert set(raw) == {0,4,-4}
        h = [v//4 for v in raw]
        D = [i for i in range(16) if h[i]]
        C = [i for i in range(16) if not h[i]]
        assert len(D)==6 and rank(16|i for i in D)==5
        assert all(rank(16|i for i in A)==5 for A in combinations(D,5))
        c_spec = walsh([int(i in C) for i in range(16)])
        assert all(abs(c)==2 for c in c_spec[1:])
        counts = Counter()
        pinned = {}
        for bits in product((-1,1),repeat=10):
            H = [0]*16
            for i, c in zip(C,bits):
                H[i]=c
            ws=walsh(H)
            small_support = tuple(i for i,c in enumerate(ws) if c)
            size = len(small_support)
            counts[size]+=1
            assert (size==16) == (sum(x<0 for x in bits)%2==0)
            combined=[x+y for x,y in zip(H,h)] + [x-y for x,y in zip(H,h)]
            assert set(combined)=={-1,1}
            global_s = [a for a,c in enumerate(walsh(combined)) if c]
            assert len(global_s)==6+size and rank(global_s)==5
            if size==7:
                tau, = [i for i,c in enumerate(ws) if abs(c)==8]
                sign=1 if ws[tau]>0 else -1
                errors=[i for i in C if H[i]!=sign*(-1)**((tau&i).bit_count()%2)]
                y, = errors
                key=(tau, small_support)
                assert pinned.setdefault(key,y)==y
                assert len([i for i in small_support if abs(ws[i])==4])==6
                assert rank([16|i for i in X]+[32|i for i in small_support])==6
                if model is None:
                    model=({i:c//4 for i,c in enumerate(walsh(h)) if c},
                           {i:c//4 for i,c in enumerate(ws) if c})
        assert counts=={7:320,10:192,16:512}
        totals.update(counts)
    return totals,model


def selector(model):
    current={1:1}
    d=1
    denominator=1
    log=[]
    for j in range(1,4):
        coeff={}
        for child, branch in enumerate(model):
            for outer, p in branch.items():
                for inner, a in current.items():
                    x=outer | (inner << (4+child*d))
                    assert x not in coeff
                    coeff[x]=a*p
        denominator*=4
        assert sum(p*p for p in coeff.values())==denominator**2
        autocorr=defaultdict(int)
        for (a,p),(b,q) in combinations(coeff.items(),2):
            autocorr[a^b]+=p*q
        assert not any(autocorr.values())
        d=rank(coeff)
        assert d==5*2**j-4 and len(coeff)==13**j and 0 not in coeff
        log.append(dict(j=j,s=len(coeff),d=d,differences=len(autocorr)))
        current=coeff
    return log


if __name__=='__main__':
    exact,model=six()
    print(json.dumps(dict(python=platform.python_version(),tolerance=0,
                          fixed_four_and_five=small(),
                          all_32_six_signs_1024_completions_each=exact,
                          repeated=selector(model)),indent=2))
