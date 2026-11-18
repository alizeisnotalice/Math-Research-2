"""Run python3 -B check_overlap_cut.py. Exhaustive exact checks for n=3,4."""
from collections import Counter
from functools import lru_cache
import json
import platform


def transform(table, n):
    values = [1-2*((table >> x)&1) for x in range(1<<n)]
    w = 1
    while w < len(values):
        for base in range(0,len(values),2*w):
            for j in range(base,base+w):
                x,y=values[j],values[j+w]
                values[j],values[j+w]=x+y,x-y
        w *= 2
    return values


@lru_cache(None)
def span(labels):
    H = frozenset({0})
    for v in labels:
        if v not in H:
            H |= frozenset(w^v for w in H)
    return H


def fourier_info(table, n):
    row = transform(table,n)
    spectrum = {v for v,x in enumerate(row) if x}
    H = span(tuple((1<<n)|v for v in sorted(spectrum)))
    return row,spectrum,H,len(H).bit_length()-1


def physical_order(g,n):
    kernel=[0]
    for x in range(1<<n):
        if (g&x).bit_count()%2==0 and x not in kernel:
            kernel += [x^z for z in kernel]
    pivot = next(x for x in range(1<<n) if (g&x).bit_count()%2)
    return kernel+[x^pivot for x in kernel]


def remap_chunks(order, n):
    picks = []
    size = 1 << n
    for offset in range(0,size,8):
        bits = min(8,size-offset)
        picks.append([sum(((val>>(x-offset))&1)<<idx for idx,x in enumerate(order)
                          if offset<=x<offset+bits) for val in range(1<<bits)])
    return picks


def audit(n):
    M = 1 << (n-1)
    tables = 1 << (2*M)
    children = [fourier_info(t,n-1) for t in range(1<<M)]
    parents = [fourier_info(t,n)[1:] for t in range(tables)]
    results=[]
    for tid in range(tables):
        w0,S0,H0,R0=children[tid % (1<<M)]
        w1,S1,H1,R1=children[tid >> M]
        Sp,_,Rp=parents[tid]
        k=len(S0 & S1)
        j=len(H0 & H1).bit_length()-1
        e=sum(abs(w0[v])==abs(w1[v]) for v in S0 & S1)
        U,L=R0+R1+2,len(S0)+len(S1)
        assert len(Sp)==2*L-2*k-e and Rp+1<=U-j
        basic=(2*L*(2*U*j-j*j) >= (2*k+e)*U*U)
        eligible=basic or k<=j
        if eligible:
            assert max(((R0+1)**2)*len(Sp)*len(S1),((R1+1)**2)*len(Sp)*len(S0)) >= (
                (Rp+1)**2*len(S0)*len(S1))
        informative = Rp > R0+R1-j
        results.append((basic,k<=j,eligible,informative))
    counts=Counter()
    any_query=set()
    for g in range(1,2*M):
        picks=remap_chunks(physical_order(g,n),n)
        seen=bytearray(tables)
        for table in range(tables):
            tid=picks[0][table & 255]
            if len(picks)==2:
                tid |= picks[1][table>>8]
            assert not seen[tid]
            seen[tid]=1
            assert len(parents[table][0])==len(parents[tid][0]) and parents[table][2]==parents[tid][2]
            basic,kj,eligible,informative=results[tid]
            counts['all_table_direction_pairs']+=1
            if informative:
                counts['informative_pairs']+=1
                if kj:
                    counts['informative_k_le_j']+=1
                if basic:
                    counts['informative_numeric_condition']+=1
                if eligible:
                    counts['informative_eligible']+=1
                    any_query.add(table)
        assert all(seen)
    counts["tables_with_informative_eligible_cut"]=len(any_query)
    return {"variables":n,"truth_tables":tables,"counts":dict(counts)}


if __name__ == "__main__":
    print(json.dumps({"python":platform.python_version(),"seed":None,"tolerance":0,
                      "results":[audit(3),audit(4)]},indent=2))
