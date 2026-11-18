"""Exact census of small affine partitions and checks of the two stated examples."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
import json
import platform

def rr(rows):
    piv = {}
    for row in rows:
        a = list(map(Fraction,row))
        for j in sorted(piv):
            c=a[j]
            if c:
                a=[x-c*y for x,y in zip(a,piv[j])]
        if any(a):
            j=next(i for i,x in enumerate(a) if x)
            piv[j]=[x/a[j] for x in a]
    return len(piv)

def br(v):
    P={}
    for x in v:
        while x.bit_length() in P:
            x ^= P[x.bit_length()]
        if x:
            P[x.bit_length()]=x
    return len(P)

def cover_analysis(n,blocks):
    N=1<<n
    idx=list(blocks)
    assert len({sum(w for s,w in blocks.items() if x in s) for x in range(N)})==1
    def norms(s):
        b=min(s)
        return [v for v in range(N) if all(not ((x^b)&v).bit_count()%2 for x in s)]
    t=br(v for s in idx for v in norms(s))
    count=max(sum(x in s for s in idx) for x in range(N))
    allrank=rr([[N*int(x in s)-len(s) for x in range(N)] for s in idx])
    local=[rr([[int(y in s) for y in range(N)] for s in idx if x not in s]) for x in range(N)]
    assert t<=min(local)<=min(len(idx)-count,allrank)
    return dict(m=len(idx),t=t,kmax=count,rank_centered=allrank,min_inactive_rank=min(local))

def overlap_example():
    P=[({1},2),({2,3},1),({0,2},1),({0,3},1)]
    R=[({x},1) for x in range(4)]+[({0,1},1),({2,3},1)]
    blocks={frozenset(x for x in range(16) if x%4 in s):w for s,w in P}
    blocks.update({frozenset(x for x in range(16) if x//4 in s):w for s,w in R})
    return cover_analysis(4,blocks)

def multiply(x,y,b,polynomial):
    value=0
    while y:
        if y&1:
            value ^= x
        y >>= 1
        x <<= 1
        if x>>b:
            x ^= polynomial
    return value

def walsh(a):
    a=list(a)
    w=1
    while w<len(a):
        for i in range(0,len(a),2*w):
            for j in range(i,i+w):
                v,t=a[j],a[j+w]
                a[j],a[j+w]=v+t,v-t
        w*=2
    return a

def family(b,polynomial):
    N=1<<b
    dataN=N*N//2
    f=[]
    fibers={}
    for u in range(N):
        row=[]
        for z in range(N//2):
            M=sum(1<<j for j in range(b) if ((multiply(1<<j,u,b,polynomial)>>1)&z).bit_count()%2)
            row.extend(1-2*(y==M) for y in range(N))
        w=walsh(row)
        assert sum(map(abs,w))*N==(3*N-4)*dataN
        for hz,hy in product(range(N//2),range(N)):
            if w[hz*N+hy]:
                fibers.setdefault(hz*N+hy,[]).append(u)
        f+=row
    S=[i for i,z in enumerate(walsh(f)) if z]
    allsets={frozenset(v) for k,v in fibers.items() if k}
    assert allsets==set(map(frozenset,combinations(range(N),2)))
    blocks={s:1 for s in allsets}
    blocks[frozenset(range(N))]=1
    cert=cover_analysis(b,blocks)
    model=dict(b=b,n=N,global_s=len(S),global_d=br(S),**cert)
    assert len(S)==1+N//2*(N*(N-1)//2)
    assert cert['min_inactive_rank']==N-1
    return model

if __name__=="__main__":
    print(json.dumps(dict(python=platform.python_version(),exact=True,
          strict_rank_example=overlap_example(),boolean_examples=[family(3,0b1011),family(4,0b10011)]),indent=2))
