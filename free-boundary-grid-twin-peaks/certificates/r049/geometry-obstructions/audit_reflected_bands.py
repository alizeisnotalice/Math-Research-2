from fractions import Fraction
from audit_reflected_seam import cluster_factors
from audit_suffix_clusters import components, connected_word, dependency, normalized_completion, square, support


def band_shape(L, factors):
    width=2*L+9
    count=4*L+1
    n=count*(width+1)-1
    c=n//8
    mid=(n-1)//2
    V,N=square(n)
    reflect=lambda v:(v[0],n-1-v[1])
    span=tuple(range(c+6,n-c-6,12))
    left=frozenset((x,y) for x,y in V
                   if x<c or (y//(width+1))%2==0 and y%(width+1)<width and x<n-c-1)
    right=frozenset((x,y) for x,y in V
                    if x>=n-c or (y//(width+1))%2==1 and y%(width+1)<width and x>c)
    assert reflect((n-c,mid))==(n-c,mid)
    assert not any(N[v]&right for v in left)
    assert support(left|right,N)==V
    holes=[set(),set()]
    for j in range(4*L):
        gap=j*(width+1)+width
        sign=-1 if j%2==0 else 1
        for x in span:
            for side, dy in ((0,sign),(1,-sign)):
                a=x+(0 if j<2*L else 2*side-1)
                cy=gap+dy
                holes[side].update(((a,cy),(a-1,cy),(a+1,cy),(a,cy+dy)))
    assert holes[0]<=left and holes[1]<=right
    A,B=left-holes[0],right-holes[1]
    Bprime=frozenset(map(reflect,B))
    w,p,q,v=(n-c-1,mid),(n-c-2,mid-1),(n-c,mid),(n-c,mid+1)
    shape_counts=[]
    states=[]
    adj2=dependency(N)
    pair_count=2*L*len(span)
    for type_i,R in enumerate((B,Bprime)):
        assert not any(N[s]&R for s in A)
        S=A|R|{w}
        assert N[w]&A and N[w]&R
        U=V-support(S,N)
        assert len(U)==8*pair_count
        obs=sorted(map(len,components(U,adj2)))
        assert obs==sorted(([3,5] if type_i==0 else [4,4])*pair_count)
        g=normalized_completion(S,V,N,adj2)
        expected=(factors[0]*factors[2] if type_i==0 else factors[1]**2)**pair_count
        assert g==expected
        aa=connected_word([p],A,N)
        bb=connected_word([q,v],R,N)
        word=connected_word([p,q,v]+aa[1:]+bb[2:]+[w],V,N)
        ranks={vertex:i for i,vertex in enumerate(word,1)}
        assert {vertex for vertex in V if all(ranks[vertex]<ranks[y] for y in N[vertex])}=={p,q}
        assert ranks[q]==2 and ranks[w]==len(S) and min(ranks[y] for y in N[q])==3
        for branch in (A,R):
            deg2={x for x in branch if len(N[x]&branch)==2}
            assert all(len(C)<=16 for C in components(deg2,N))
        shape_counts.append(g)
        states.append((len(A),len(R),len(U),len(support(S,N)-S)))
    assert states[0]==states[1]
    assert shape_counts[0]/shape_counts[1]==Fraction(20300,18513)**pair_count
    return n,width,pair_count,states[0]


if __name__=="__main__":
    factors=cluster_factors()
    for L in (1,2,3,5):
        print(f"L={L}: n, band width, paired motifs, shape counts =",band_shape(L,factors))
