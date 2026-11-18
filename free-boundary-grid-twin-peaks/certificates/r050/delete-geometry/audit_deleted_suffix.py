from collections import defaultdict
from fractions import Fraction


def neighbors_for(width, height):
    points=[(x,y) for y in range(height) for x in range(width)]
    lookup={p:i for i,p in enumerate(points)}
    neighbors=[]
    for x,y in points:
        neighbors.append(sum(1<<lookup[v] for v in ((x-1,y),(x+1,y),(x,y-1),(x,y+1)) if v in lookup))
    return points, lookup, neighbors


def take(mask):
    while mask:
        low=mask & -mask
        yield low.bit_length()-1
        mask ^= low


def union_sides(mask, neighbors):
    value=0
    for v in take(mask):
        value |= neighbors[v]
    return value


def is_connected(mask,neighbors):
    if not mask:
        return False
    state=mask & -mask
    frontier=state
    while frontier:
        frontier=union_sides(frontier,neighbors) & mask & ~state
        state |= frontier
    return state==mask


def completion_counts(neighbors, blocked=0):
    N=len(neighbors)
    full=((1<<N)-1)^blocked
    h=[0]*(1<<N)
    boundary=[0]*(1<<N)
    h[full]=1
    for mask in range(1,full+1):
        low=mask & -mask
        boundary[mask]=boundary[mask^low]|neighbors[low.bit_length()-1]
    for S in range(full-1,0,-1):
        if S&blocked:
            continue
        candidates=boundary[S] & (full^S)
        h[S]=sum(h[S|(1<<v)] for v in take(candidates))
    return h,boundary


def census(n):
    pts,_,adj=neighbors_for(n,n)
    H,boundary=completion_counts(adj)
    full=(1<<len(pts))-1
    records=defaultdict(lambda:[0,Fraction(1),None,0])
    for b in (1,2):
        for B in range(1,full+1):
            if B.bit_count()!=b or not is_connected(B,adj):
                continue
            assert is_connected(full^B,adj)
            D,_=completion_counts(adj,B)
            allowed=full & ~(B|boundary[B])
            for w in take(boundary[B]&~B):
                subset=allowed
                while subset:
                    if (adj[w]&subset) and is_connected(subset,adj):
                        S=subset|(1<<w)|B
                        K=S^B
                        assert H[S]>0 and D[K]>0
                        ratio=Fraction(D[K],H[S])
                        bag=records[b]
                        bag[0]+=1
                        if ratio<bag[1]:
                            bag[1]=ratio
                            bag[2]=(tuple(pts[v] for v in take(B)),pts[w],tuple(pts[v] for v in take(subset)))
                        R=full &~S
                        U=R&~boundary[S]
                        extras=R&boundary[B]&~boundary[K]
                        protected=R&boundary[K]
                        # Check the bound under the staged graphical condition.
                        test=True
                        old=U
                        for d in take(extras):
                            options=[(z,t) for z in take(adj[d]&protected) for t in take(protected)
                                     if z!=t and not (adj[d]&old&~adj[t])]
                            if not options:
                                test=False
                                break
                            old|=1<<d
                        if test:
                            assert ratio>=Fraction(1,2**extras.bit_count())
                            bag[3]+=1
                    subset=(subset-1)&allowed
    for b,record in records.items():
        total,value,shape,gate=record
        print(f"n={n}, b={b}: {total} legal triples, min ratio={value}, witness={shape}, detour rule qualified={gate}")


def strip(n):
    # The residual ladder has n-1 rows, x=0 one extra vertex for deleting B.
    pts,_,adj=neighbors_for(2,n-1)
    N=len(pts)
    full=(1<<N)-1
    allowed=sum(1<<i for i,v in enumerate(pts) if v[0]==1)
    def words(free):
        dp=[0]*(full+1)
        dp[0]=1
        for mask,value in enumerate(dp):
            if value:
                for v in take(full^mask):
                    if ((1<<v)&free) or (adj[v]&mask):
                        dp[mask|(1<<v)]+=value
        return dp[full]
    p=Fraction(words(allowed),words(allowed|1))
    assert p>=Fraction(1,2)
    return p


def threshold_gate():
    tested=0
    for ai in range(21):
        for bi in range(21):
            a,b=Fraction(ai,21),Fraction(bi,21)
            A=1-max(a,b)
            B=max(b-a,0)*(1-b)
            after=A*(b+(1-b)**2/2)+B*(1-b)/2
            assert 2*after>=A+B
            tested+=1
    print("exact threshold gate checks:",tested)


def zero_family():
    for n in range(4,13):
        pts,ids,adj=neighbors_for(n,n)
        full=(1<<len(pts))-1
        B=sum(1<<ids[v] for v in ((1,0),(1,1),(0,1)))
        w=ids[(2,1)]
        for A in ({(3,1)}, {v for v in pts if v[0]>=3}):
            As=sum(1<<ids[v] for v in A)
            K=As|(1<<w)
            S=K|B
            assert is_connected(As,adj) and is_connected(B,adj)
            assert not (union_sides(B,adj)&As)
            assert adj[w]&B and adj[w]&As and is_connected(S,adj)
            assert not is_connected(full^B,adj)
            corner=ids[(0,0)]
            assert not (adj[corner]&(full^B))
    print("zero probability construction: 18 thin/bulk triples")


if __name__=="__main__":
    for n in (3,4):
        census(n)
    for n in range(3,10):
        print(f"n={n}, integrated two-column strip ratio={strip(n)}")
    threshold_gate()
    zero_family()
