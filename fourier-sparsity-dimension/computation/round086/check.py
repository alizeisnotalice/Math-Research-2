"""Exact finite audits for the external mass gap and seed substitution."""
from collections import Counter
from itertools import product
from random import Random
import json
import platform


def walsh(data):
    arr = list(data)
    step = 1
    while step < len(arr):
        for i in range(0, len(arr), 2*step):
            for j in range(i, i+step):
                x,y = arr[j],arr[j+step]
                arr[j],arr[j+step] = x+y,x-y
        step *= 2
    return arr


def rank(vectors):
    piv = {}
    for v in vectors:
        while v:
            i = v.bit_length()
            if i not in piv:
                piv[i]=v
                break
            v ^= piv[i]
    return len(piv)


def tag(row):
    dark = [x for x,val in enumerate(row) if val < 0]
    if len(dark) in (0,len(row)):
        return 'constant'
    if len(dark) in (1,len(row)-1):
        return 'single'
    return 'other'


def one_core_bit():
    # By a translation and sign, a full four-dimensional active slice starts here.
    init = [-1]+[1]*15
    a = walsh(init)
    tally = Counter()
    masses = Counter()
    for mask in range(1<<16):
        row = [1-2*((mask>>i)&1) for i in range(16)]
        b = walsh(row)
        # Two-point Walsh triangle identity gives exact global Fourier mass times 16.
        outside = sum(max(abs(x),abs(y)) for x,y in zip(a[1:],b[1:]))
        if outside < 40:
            kind = tag(row)
            assert kind in ('constant','single') and outside == 30
            tally[kind] += 1
        masses[outside] += 1
    return dict(tested=sum(masses.values()),accepted=dict(tally),
                next_mass_numerator_over_16=min(v for v in masses if v>30))


def on_two_core_bits():
    # Each row is constant (+ or -) or an active singleton at one of 16 locations.
    options = [(0,None),(1,None)] + list(product(range(2),range(16)))
    tally = Counter()
    smallest_other = 10**9
    for tail in product(options,repeat=3):
        cfg = [(0,0),*tail]
        active = [y for y,(_,x) in enumerate(cfg) if x is not None]
        affine = len(active) in (1,2,4)
        aligned = len(active)<4 or (
            sum(s for s,_ in cfg)%2 == 0 and cfg[0][1]^cfg[1][1]^cfg[2][1]^cfg[3][1] == 0)
        good = affine and aligned
        # For each nonzero label, core coefficient is +-1/8 times ternary indicator.
        total = 0
        for beta in range(1,16):
            col = [(1-2*(s^((x&beta).bit_count()%2))) if x is not None else 0 for s,x in cfg]
            total += sum(map(abs,walsh(col)))
        if total < 80:
            assert good and total == 60
            tally['affine'] += 1
        else:
            assert not good
            tally['nonaffine'] += 1
            smallest_other = min(smallest_other,total)
    return dict(tested=sum(tally.values()),counts=dict(tally),next_mass_numerator_over_32=smallest_other)


def variable_tree(x):
    # Three-leaf binary tree. First two coordinates are internal nodes.
    leaf = 0 if not x&1 else 1+((x>>1)&1)
    return 1-2*((x>>(2+leaf))&1)


def compare_mutations():
    rng = Random(86041)
    seeds = []
    for k in (4,5):
        size = 1<<k
        if k==4:
            size2=32
            j=[variable_tree(x) for x in range(32)]
        else:
            size2=128
            j=[]
            for x in range(128):
                leaf = 2*((x&1)!=0) + ((x>>(1+int(bool(x&1))))&1)
                j.append(1-2*((x>>(3+leaf))&1))
        anchor=next(i for i,c in enumerate(walsh(j)) if c)
        j=[v*(1-2*((x&anchor).bit_count()%2)) for x,v in enumerate(j)]
        sj=sum(c!=0 for c in walsh(j))
        assert sj==(10 if k==4 else 16)
        for _ in range(24):
            normal=rng.randrange(1,8)
            rhs=rng.randrange(2)
            active=[(y&normal).bit_count()%2==rhs for y in range(8)]
            phase=rng.randrange(8)
            outside=[rng.choice([-1,1]) for _ in range(8)]
            old=[]; new=[]
            for y,on in enumerate(active):
                sign=1-2*((y&phase).bit_count()%2)
                old.extend(sign*(1-2*int(x==0)) if on else outside[y] for x in range(size))
                new.extend(sign*j[x] if on else outside[y] for x in range(size2))
            oldspec=[i for i,c in enumerate(walsh(old)) if c]
            newspec=[i for i,c in enumerate(walsh(new)) if c]
            ro=rank(len(old)|v for v in oldspec)
            rn=rank(len(new)|v for v in newspec)
            assert len(newspec)<=len(oldspec)-((size-1)-sj)*2 and rn>=ro
            seeds.append((k,len(oldspec),ro,len(newspec),rn))
    return dict(tested=len(seeds),examples=sorted(set(seeds)))


if __name__=='__main__':
    print(json.dumps(dict(python=platform.python_version(),tolerance=0,seed=86041,
                          core_one=one_core_bit(),core_two=on_two_core_bits(),
                          substitutions=compare_mutations()),indent=2))
