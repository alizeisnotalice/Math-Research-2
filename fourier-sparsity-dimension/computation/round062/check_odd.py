"""Exact small Boolean phase-retrieval and extremal denominator checks."""
import json
import platform
from collections import defaultdict


def walsh(values):
    row=list(values)
    size=1
    while size<len(row):
        for i in range(0,len(row),2*size):
            for j in range(i,i+size):
                x,y=row[j],row[j+size]
                row[j],row[j+size]=x+y,x-y
        size*=2
    return row


def rank(vectors):
    basis={}
    for v in vectors:
        while v:
            i=v.bit_length()
            if i not in basis:
                basis[i]=v
                break
            v ^= basis[i]
    return len(basis)


def truth(mask,N):
    return [1-2*((mask>>i)&1) for i in range(N)]


def translate_mask(mask,N,b):
    return sum(((mask>>(x^b))&1)<<x for x in range(N))


def census(a):
    N=1<<a
    classes=defaultdict(set)
    extremal=0
    for mask in range(1<<N):
        w=walsh(truth(mask,N))
        support=[i for i,z in enumerate(w) if z]
        div=min(abs(w[i])&-abs(w[i]) for i in support)
        k=a-(div.bit_length()-1)
        if k:
            odd=[i for i in support if (w[i]//div)%2]
            assert len(odd)>=2**(k+1)
            if len(support)==2**(k+1):
                extremal+=1
                assert set(odd)==set(support)
                assert rank(x^support[0] for x in support)==k+1
        if mask.bit_count()%2:
            classes[tuple(map(abs,w))].add(mask)
    for members in classes.values():
        seed=next(iter(members))
        orbit={translate_mask(seed,N,b)^c for b in range(N) for c in (0,(1<<N)-1)}
        assert orbit==members
    return dict(a=a,tables=1<<N,odd_tables=sum(map(len,classes.values())),
                magnitude_classes=len(classes),extremal_granularity_tables=extremal)


def padding(a,seed):
    original=truth(seed,8)
    base=walsh(original)
    sg=sum(v!=0 for v in base)
    data=[]
    for choose in range(2):
        for y in range(1<<a):
            for v in range(2):
                for x in range(8):
                    data.append((1-2*(y==0)) if not choose else (1-2*v)*original[x])
    w=walsh(data)
    W=1<<(a+1)
    support=[i for i,y in enumerate(w) if y]
    unit=len(w)//(1<<a)
    O=[i for i in support if w[i]//unit%2]
    assert all(w[i]%unit==0 for i in support)
    assert set(O)=={i for i in range(len(w)) if not i%16}
    assert len(support)==W+2*sg and len(O)==W
    assert rank(support)>=rank(i for i,x in enumerate(base) if x)
    return dict(a=a,seed=seed,s_G=sg,M=W,s_F=len(support),d_F=rank(support))


print(json.dumps(dict(python=platform.python_version(),
    small=[census(n) for n in (2,3,4)],
    padded=[padding(a,seed) for a in (4,5,7) for seed in (0x16,0x69,0x81)]),indent=2))
