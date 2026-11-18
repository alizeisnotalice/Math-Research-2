from collections import deque
from fractions import Fraction
from itertools import product
from math import factorial


def graph(n):
    V=frozenset(product(range(n),repeat=2))
    L={v:V&{(v[0]-1,v[1]),(v[0]+1,v[1]),(v[0],v[1]-1),(v[0],v[1]+1)} for v in V}
    return V,L


def flood(seeds, vertices, N):
    q=deque(seeds)
    connected=set(seeds)
    while q:
        v=q.popleft()
        new=(N[v]&vertices)-connected
        q.extend(new)
        connected |= new
    return connected


def shape(n,A,B,w):
    V,N=graph(n)
    A,B=set(A),set(B)
    K=A|{w}
    S=K|B
    assert flood({next(iter(A))},A,N)==A
    assert flood({next(iter(B))},B,N)==B
    assert not any(N[x]&B for x in A) and N[w]&A and N[w]&B
    H=V-S
    free={x for x in H if N[x]&K}
    free_orig={x for x in H if N[x]&S}
    return H,N,K,free,free_orig


def completion(H,N,free):
    points=sorted(H)
    ids={v:i for i,v in enumerate(points)}
    avail=[v in free for v in points]
    masks=[sum(1<<ids[w] for w in N[v]&H) for v in points]
    full=(1<<len(points))-1
    dp=[0]*(full+1)
    dp[0]=1
    for s,value in enumerate(dp):
        if value:
            unused=full^s
            while unused:
                bit=unused&-unused
                i=bit.bit_length()-1
                unused^=bit
                if avail[i] or masks[i]&s:
                    dp[s|bit]+=value
    return dp[full]


def companion_case(n,b,exact=False):
    A={(3,y) for y in range(1,b+2)}|{(2,b+1)}
    B={(1,y) for y in range(1,b+1)}
    w=(2,1)
    H,N,K,F,F_orig=shape(n,A,B,w)
    U=H-F_orig
    D=F_orig-F
    route=[(2,0),(1,0)]+[(0,y) for y in range(b+1)]
    P=set(route)
    d,t=(0,b),(1,b+1)
    assert D<=P<=H and d in D and not (N[d]&F)
    assert all((N[v]&(U-P))==({(0,b+1)} if v==d else set()) for v in P)
    assert t in F-P and N[t]&{(0,b+1)}
    seen=set(K)
    for v in route:
        assert N[v]&seen
        seen.add(v)
    if exact:
        h1=completion(H,N,F)
        h0=completion(H,N,F_orig)
        ratio=Fraction(h1,h0)
        assert ratio>=Fraction(1,2*factorial(b+3))
        print("large-core grouped shape, n b exact target/source:",n,b,h1,h0,ratio)
    return len(U)


def embedded_case(n,b):
    ext={(x,y) for x in range(n) for y in range(n) if x>=4 or y>=4}
    if b==1:
        B={(2,0)}
        w=(2,1)
        A={(3,1),(3,2),(2,3),(3,3)}
        extra=(2,2)
        ratios=(Fraction(266,491),Fraction(47,95))
    else:
        B={(1,0),(2,0)}
        w=(3,0)
        A={(3,1),(1,3),(3,3)}
        extra=(3,2)
        ratios=(Fraction(382,2275),Fraction(11,70))
    entries=[]
    for addon,expected in zip((set(),{extra}),ratios):
        H,N,_,F,F0=shape(n,A|ext|addon,B,w)
        h1,h0=completion(H,N,F),completion(H,N,F0)
        assert Fraction(h1,h0)==expected
        entries.append((h1,h0))
    assert ratios[1]<ratios[0]
    return entries


def group_moment(k, sizes):
    # k-integral of x^(k-1) product(1-x^r).
    terms={0:1}
    for r in sizes:
        new=terms.copy()
        for degree,coefficient in terms.items():
            new[degree+r]=new.get(degree+r,0)-coefficient
        terms=new
    return sum(Fraction(k*value,k+degree) for degree,value in terms.items())


def cap_checks():
    checked=0
    for size_count in range(1,5):
        for sizes in product(range(2,7),repeat=size_count):
            k=sum(r-1 for r in sizes)
            moment=group_moment(k,sizes)
            cutoff=Fraction(1,2**k)
            for r in sizes:
                cutoff*=1-Fraction(1,2**r)
            assert 0<cutoff<=moment<=1
            checked+=1
    return checked


if __name__=="__main__":
    for n in (4,5,8,12,40):
        for b in (1,2):
            companion_case(n,b,exact=(n<=5))
    for n in (5,6,12,32):
        for b in (1,2):
            print("embedded n,b old/add-one pairs (target,source):",n,b,embedded_case(n,b))
    print("exact max-cap moment comparisons:",cap_checks())
