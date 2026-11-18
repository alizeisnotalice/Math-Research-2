"""Exact data-basis marginals on genuinely overlapping six-circuit completions."""
from fractions import Fraction as Q
from itertools import combinations
from math import comb
from random import Random
import json

def insert(basis, x):
    while x and x.bit_length() in basis:
        x ^= basis[x.bit_length()]
    if x:
        basis[x.bit_length()] = x
        return True
    return False

def fwht(v):
    v=list(v)
    step=1
    while step<len(v):
        for i in range(0,len(v),2*step):
            for j in range(i,i+step):
                x,y=v[j],v[j+step]
                v[j],v[j+step]=x+y,x-y
        step*=2
    return v

def run(seed, identical):
    rng=Random(seed); bits=5; size=1<<bits
    addresses=[i for i in range(16) if (i%4==0)==(i//4==0)]
    slices=[]; choices=[]
    for i in range(10):
        if identical or i<4: support=[1,2,4,7]
        elif i<8:
            while True:
                v=rng.sample(range(1,16),3)
                if v[0]^v[1]!=v[2]: break
            support=v+[v[0]^v[1]^v[2]]
        else: support=[1 if i==8 else 16]
        values={b:1 for b in support}
        if len(support)==4: values[rng.choice(support)]=-1
        else: values[support[0]]=2
        slices.append(values)
        choices.append(list(combinations(support,3)) if len(support)==4 else [tuple(support)])
    tables=[]
    for g in slices:
        v=fwht([g.get(x,0) for x in range(size)])
        assert set(v)<= {2,-2}
        tables.append([x//2 for x in v])
    truth=[]
    for x in range(size):
        for t in range(2):
            for y in range(16):
                h=int(y%4==0)-int(y//4==0)
                truth.append(h*(-1)**t if h else tables[addresses.index(y)][x])
    spectrum=fwht(truth)
    lifted={b:sum(bool(z) for z in spectrum[b*32:b*32+32]) for b in range(1,size)}
    before=[]
    for mask in range(1024):
        basis={}
        for i,g in enumerate(slices):
            if mask>>i&1:
                for x in g: insert(basis,x)
        before.append(basis)
    ps=[Q(0)]*size; average_bound=[Q(0)]*size
    source_cost=Q(0)
    for i, options in enumerate(choices):
        marginal=Q(len(options[0]),len(slices[i]))
        source_cost+=len(slices[i])*marginal**2/16
        for b in slices[i]:
            average_bound[b]+=marginal/sum(b in g for g in slices)
        for mask,basis in enumerate(before):
            if mask>>i&1: continue
            mass=Q(1,10*comb(9,mask.bit_count())*len(options))
            for candidate in options:
                pivots=basis.copy()
                for b in candidate:
                    if insert(pivots,b): ps[b]+=mass
    assert all(x<=y for x,y in zip(ps,average_bound))
    assert sum(ps)==len(before[-1])
    cost=sum(ps[b]**2/m for b,m in lifted.items() if m)
    assert cost<=source_cost
    return dict(seed=seed,identical_supports=identical,parent_s=sum(bool(x) for x in spectrum),
                projection_rank=len(before[-1]),data_cost=str(cost),budget=str(source_cost),
                max_label_multiplicity=max(sum(b in g for g in slices) for b in range(1,size)))
print(json.dumps([run(79,True)]+[run(i,False) for i in range(79,83)],indent=2))
