from collections import deque
from fractions import Fraction as F
from itertools import permutations, product
from math import factorial


def points(n):
    V=set(product(range(n),repeat=2))
    N={v:V&{(v[0]-1,v[1]),(v[0]+1,v[1]),(v[0],v[1]-1),(v[0],v[1]+1)} for v in V}
    return V,N


def build(V,N,A,B,w):
    K=A|{w}
    S=K|B
    H=V-S
    F={v for v in H if N[v]&K}
    V0={v for v in H if N[v]&S}
    U=H-V0
    D=V0-F
    assert not any(N[v]&B for v in A) and N[w]&A and N[w]&B
    return H,F,U,D


def finish(prefix, vertices, N, source):
    q=deque(prefix+sorted(source-set(prefix)))
    used=set(q)
    out=list(q)
    while q:
        v=q.popleft()
        nxt=sorted((N[v]&vertices)-used)
        q.extend(nxt)
        out.extend(nxt)
        used.update(nxt)
    assert used==vertices
    return out


def exact_old_table(H, N, source):
    seq=sorted(H)
    ids={v:i for i,v in enumerate(seq)}
    masks=[sum(1<<ids[u] for u in N[v]&H) for v in seq]
    free=[v in source for v in seq]
    full=(1<<len(H))-1
    pre=[0]*(full+1)
    suf=[0]*(full+1)
    pre[0]=suf[full]=1
    def edges(S):
        left=full^S
        while left:
            bit=left&-left
            left^=bit
            i=bit.bit_length()-1
            if free[i] or masks[i]&S:
                yield bit
    for s in range(full+1):
        if pre[s]:
            for edge in edges(s):
                pre[s|edge]+=pre[s]
    for s in range(full-1,-1,-1):
        if pre[s]:
            suf[s]=sum(suf[s|edge] for edge in edges(s))
    assert suf[0]==pre[full]
    return ids,pre,suf,masks,free


def audit_witness_pool():
    # Opposite witnesses u,v have separate safe backups t,q; z is an extra protected neighbor.
    d,z,t,q,u,v,x,y=range(8)
    pair={(d,z),(d,u),(d,v),(t,u),(q,v),(x,u),(y,v),(x,z),(y,z)}
    N={r:{b if a==r else a for a,b in pair if r in (a,b)} for r in range(8)}
    C={u,v,x,y}
    before=after=s0_sum=best_sum=rank2_sum=0
    for word in permutations(range(8)):
        R={p:j for j,p in enumerate(word,1)}
        def peak(k):
            return all(R[k]<R[n] for n in N[k])
        if any(map(peak,C)):
            continue
        before+=1
        after+=not peak(d)
        W={u0 for u0 in N[d]&C if all(R[w]>R[u0] for w in N[u0]-{d})}
        js=[min((R[u0] for u0 in W-N[t0]),default=9) for t0 in (t,q)]
        ranked=sorted(R[u0] for u0 in W)
        j2=ranked[1] if len(ranked)>1 else 9
        assert max(js)>=j2
        W0={u0 for u0 in N[d]&C if all(R[w]>R[u0] for w in N[u0]-{d,t,q,z})}
        pick=max((t,q),key=lambda t0:min((R[u0] for u0 in W0-N[t0]),default=9))
        adaptive=min((R[u0] for u0 in W-N[pick]),default=9)
        s0_sum+=adaptive
        best_sum+=max(js)
        rank2_sum+=j2
    assert 18*after>=s0_sum
    assert 36*after>=best_sum>=rank2_sum
    return before,after,s0_sum,best_sum,rank2_sum


def audit_grid_pool():
    for n in (7,8,12,30):
        V,N=points(n)
        shift=lambda p:(p[0]+2,p[1]+1)
        A={(-2,y) for y in range(-1,4)}|{(x,3) for x in range(-2,4)}|{(3,y) for y in range(-1,4)}|{(2,1),(2,-1)}
        H,F,U,D=build(V,N,{shift(v) for v in A},{shift((1,0))},shift((2,0)))
        d,z,u,v,t,q=map(shift,((0,0),(-1,0),(0,1),(0,-1),(-1,1),(-1,-1)))
        assert D=={d} and N[d]&U=={u,v} and {z,t,q}<=F and z in N[d]
        assert N[u]&N[v]=={d}
        for prefix,expected in (([d,u,q,v],{u}),([d,u,v],{u,v})):
            word=finish(prefix,H,N,F|D)
            R={a:i for i,a in enumerate(word,1)}
            out={a for a in N[d]&U if {n for n in N[a]&H if R[n]<R[a]}=={d}}
            assert out==expected and R[u]==2
    return 8


def audit_prefix_shape():
    V,N=points(4)
    H,protected,U,D=build(V,N,{(2,2)},{(1,1)},(2,1))
    P=((2,0),(1,0),(0,0),(0,1))
    ps=set(P)
    assert D=={(1,0),(0,1)} and D<=ps<=H
    assert (2,0) in protected
    assert all(a in N[b] for a,b in zip(P,P[1:]))
    surrogate=(1,2)
    assert surrogate in protected and N[P[-1]]&(U-ps)=={(0,2)}
    assert (0,2) in N[surrogate]
    outer={v for p in P[:-1] for v in N[p]&(U-ps)}
    assert outer=={(3,0)}
    u=(3,0)
    ids,pre,suf,masks,free=exact_old_table(H,N,protected|D)
    total=pre[-1]
    after=exact_old_table(H,N,protected)[1][-1]
    flag=1<<ids[u]
    frozen=sum(1<<ids[x] for x in (N[u]&H)-(ps|{surrogate}))
    count=rank=0
    for S,C in enumerate(pre):
        if not C or S&(flag|frozen):
            continue
        if masks[ids[u]]&S or free[ids[u]]:
            v=C*suf[S|flag]
            count+=v
            rank+=v*(S.bit_count()+1)
    assert 0<count<total
    m=len(H)
    bound=F(total-count,1536)+F(rank,768*(m+1))
    assert after>=bound
    return m,total,after,count,rank


def audit_cap_constants():
    checked=0
    for k in range(1,7):
        for sizes in ((),(2,),(2,3),(3,4,4)):
            g=len(sizes)
            for q in (F(0),F(1,9),F(1,2),F(2,3),F(9,10)):
                z=max(q,F(1,2))
                for j in range(k+1):
                    for c in range(g+1):
                        upper=(1+z)/2
                        exact=F(1,factorial(j)*factorial(k-j))*(upper-z)**j/2**(k-j)
                        for t,r in enumerate(sizes):
                            exact *= 1-(upper if t<c else F(1,2))**r
                        lower=(1-q)**(j+c)*F(1,2**(k+j+g+c)*factorial(j)*factorial(k-j))
                        assert exact>=lower
                        checked+=1
    return checked


if __name__=="__main__":
    print("8-point whole-integral witness pool counts:",audit_witness_pool())
    print("grid unique/both early witness words:",audit_grid_pool())
    print("actual 4x4 prefix repair (m,old,del,active,active rank sum):",audit_prefix_shape())
    print("rational prefix-volume comparisons:",audit_cap_constants())
