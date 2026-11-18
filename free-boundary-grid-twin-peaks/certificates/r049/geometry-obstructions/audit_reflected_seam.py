from fractions import Fraction
from itertools import combinations
from audit_suffix_clusters import components, connected_word, dependency, local_nopeaks, normalized_completion, square, support


def cluster_factors():
    V,N=square(13)
    records=[]
    for delta in (0,1,2):
        C={(5,5),(7,5+delta)}|{(6,y) for y in range(4,7) if abs(y-(5+delta))<=1}
        local=local_nopeaks(C,N)
        J={frozenset():Fraction(1)}
        sums=[Fraction(1)]
        for k in range(1,len(C)+1):
            summed=Fraction(0)
            for I0 in combinations(sorted(C),k):
                I=frozenset(I0)
                if any(N[v]&I for v in I):
                    continue
                J[I]=sum(J[I-{v}] for v in I)/len(support(I,N))
                summed += J[I]
            sums.append(summed)
        assert local == sum((-1)**k*term for k,term in enumerate(sums))
        records.append(local)
        print(f"delta={delta}: no peaks={local}, inclusion-exclusion layer sums={sums}")
    assert records[0]*records[2]/records[1]**2==Fraction(20300,18513)
    return records


def seam_shape(K):
    n,m=24*K+13,12*K+6
    V,N=square(n)
    ref=lambda v:(v[0],n-1-v[1])
    left_rows=[]
    right_rows=[]
    for j in range(K):
        row=6+12*j
        left_rows += [row,n-row]
        right_rows += [row,n-2-row]
    def holes_at(x, rows, dx):
        return {(x,y+s) for y in rows for s in (-1,0,1)}|{(x+dx,y) for y in rows}
    A=frozenset(v for v in V-holes_at(m-1,left_rows,-1) if v[0]<m)
    B=frozenset(v for v in V-holes_at(m+1,right_rows,1) if v[0]>m)
    Bright=frozenset(map(ref,B))
    w,p,q=(m,m),(m-1,m-1),(m+1,m)
    assert ref(w)==w and ref(q)==q
    assert len(A)==len(B)==len(Bright)==m*n-8*K
    assert all({ref(x) for x in N[v]&B} == N[ref(v)]&Bright for v in B)
    assert len(components(A,N))==len(components(B,N))==len(components(Bright,N))==1
    test=[]
    for side in (B,Bright):
        assert not any(N[a] & side for a in A)
        S=A|side|{w}
        assert len(components(S,N))==1 and len(components(S-{w},N))==2
        wordA=connected_word([p],A,N)
        r=(m+1,m+1)
        wordB=connected_word([q,r],side,N)
        word=connected_word([p,q,r]+wordA[1:]+wordB[2:]+[w],V,N)
        rank={v:i for i,v in enumerate(word,1)}
        assert {v for v in V if all(rank[v]<rank[u] for u in N[v])}=={p,q}
        assert rank[w]==len(S) and min(rank[u] for u in N[q])==3 and rank[q]==2
        U=V-support(S,N)
        assert len(U)==8*K
        f=normalized_completion(S,V,N,dependency(N))
        test.append(f)
    expected=Fraction(20300,18513)**K
    assert test[0]/test[1]==expected
    return n,K,len(U),4*K,len(A)


if __name__=="__main__":
    cluster_factors()
    print("shape audit (n, paired locations, undominated, per-side holes/2, side area):",
          [seam_shape(K) for K in (1,2,5,10)])
