"""Exact verification of both extremal local costs, Boolean completions and every saturated circuit pattern."""
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations
import json,platform

def walsh(a):
    a=list(a);j=1
    while j<len(a):
        for i in range(0,len(a),2*j):
            for x in range(i,i+j):
                a[x],a[x+j]=a[x]+a[x+j],a[x]-a[x+j]
        j*=2
    return a

def splits(indices):
    if not indices:
        yield []
        return
    head=indices[0]
    for size in range(2,len(indices)+1,2):
        for tail in combinations(indices[1:],size-1):
            group=(head,)+tail
            rest=[i for i in indices if i not in group]
            for ss in splits(rest):
                yield [group]+ss

def validate(K):
    p=K+1;size=1<<p
    frac=[]
    for blocks in splits(list(range(K))):
        extra=[sum(1<<i for i in block) for block in blocks]
        cols=[1<<i for i in range(K)]+extra
        assert max(sum((q&v).bit_count()%2 for v in cols) for q in range(1<<K))==K
        frac.append(sum(Q(len(B)**2,len(B)+1) for B in blocks))
    assert max(frac)==Q(K*K,K+1)
    vals=[Q(p-2*v.bit_count()-(-1)**v.bit_count(),K) for v in range(size)]
    alphabet={0,1,-1} if K==4 else {0,Q(1,2),-Q(1,2),1,-1}
    assert set(vals)==alphabet
    zn=1 if K==4 else 5
    T=[]
    private_cost=[]
    support=set()
    for g in vals:
        if K==4:
            row=[int(g+((-1)**z)*(1-g*g)) for z in range(2)]
            c=Q(int(not g))
        else:
            row=[]
            for z in range(32):
                x1,x2,x3,x4,t=((z>>i)&1 for i in range(5))
                if abs(g)==1: val=int(g)
                elif abs(g)==Q(1,2): val=(1 if g>0 else -1)*(-1)**(x1*x2)
                else: val=(-1)**(t+x3*x4)
                row.append(val)
            c={0:Q(9,4),Q(1,2):Q(4,3),1:Q(0)}[abs(g)]
        assert set(row)<={1,-1} and Q(sum(row),len(row))==g
        support |= {i for i,z in enumerate(walsh(row)) if z}
        T.extend(row)
        private_cost.append(c)
    W=walsh(T);q0=[Q(z,len(T)) for z in W[::1<<zn]]
    assert {i for i,z in enumerate(q0) if z}=={1<<j for j in range(p)}|{size-1}
    assert q0[-1]==Q(-1,K)
    pairs=Counter(i^j for i,j in combinations(support,2))
    assert all(min(pairs[a^b],pairs[a^c],pairs[a^d])==2
               for a,b,c,d in combinations(support,4) if a^b^c^d==0)
    v=max(frac)+sum(private_cost)/size
    assert v==(Q(153,40) if K==4 else Q(40751,4608))
    return dict(K=K,rank_increment=K,relative_basis_cost=str(max(frac)),
     private_cost=str(sum(private_cost)/size),joint=str(v),full_truth_table_entries=len(T),
     label_set_in_private_quotient=sorted(support),saturated_even_set_partitions=len(frac),
     zero_band_value_histogram={str(k):v for k,v in sorted(Counter(vals).items())})
print(json.dumps(dict(python=platform.python_version(),cases=[validate(4),validate(8)]),indent=2))
