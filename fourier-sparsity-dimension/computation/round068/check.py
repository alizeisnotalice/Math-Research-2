"""Small genuine Boolean clique-stars, checked by full truth tables and exact Walsh transforms."""
from itertools import combinations
from collections import Counter
import json,platform

def transform(a):
    a=list(a);k=1
    while k<len(a):
        for i in range(0,len(a),2*k):
            for j in range(i,i+k):
                a[j],a[j+k]=a[j]+a[j+k],a[j]-a[j+k]
        k*=2
    return a

def rank(X):
    p={}
    for x in X:
        while x.bit_length() in p:
            x^=p[x.bit_length()]
        if x:p[x.bit_length()]=x
    return len(p)

def case(h,ms):
    n=sum(ms);N=n//2+1;L=1<<(h+n)
    assert N&(N-1)==0
    pairs=[];owner={};off=0
    for i,m in enumerate(ms):
        u=i+1
        for j in range(m):
            e=1<<(h+off+j)
            owner[e]=owner[e^u]=i
        for j in range(0,m,2):
            pairs.append((u,1<<(h+off+j),1<<(h+off+j+1)))
        off+=m
    table=[2*int(x%(1<<h)==0)-1 for x in range(L)]
    for u,e0,e1 in pairs:
        table.extend(1-2*int(bool(x&(e1 if (x&u).bit_count()%2 else e0))) for x in range(L))
    w=transform(table); S=[i for i,x in enumerate(w) if x]
    B={i%L for i in S};C=set(range(1<<h))
    assert B==C|set(owner)
    multiplicities=Counter(a^b for a,b in combinations(B,2))
    for a,b in combinations(B,2):
        assert (multiplicities[a^b]>1)==(a in C or b in C or owner[a]==owner[b])
    assert all(sum(i%L==b for i in S)==N for b in B)
    q=sum(i%L not in C for i in S)
    assert q==n*(n+2) and rank(B)==h+n
    return dict(h=h,private_ranks=ms,truth_table_size=len(w),spectral_sparsity=len(S),
                projected_rank=rank(B),nonzero_quotient_tail=q,clique_star_pairs_checked=len(B)*(len(B)-1)//2)
print(json.dumps(dict(python=platform.python_version(),cases=[case(2,[2,4]),case(3,[2,2,2])]),indent=2))
