from collections import deque
from fractions import Fraction
from math import factorial
from random import Random


def square(n):
    vertices = frozenset((x,y) for x in range(n) for y in range(n))
    neighbors = {
        (x,y): vertices & {(x-1,y),(x+1,y),(x,y-1),(x,y+1)}
        for x,y in vertices
    }
    return vertices, neighbors


def components(points, adjacency):
    left = set(points)
    answer = []
    while left:
        start = left.pop()
        queue = deque([start])
        component = {start}
        while queue:
            v = queue.popleft()
            new = adjacency[v] & left
            left -= new
            component |= new
            queue.extend(new)
        answer.append(frozenset(component))
    return answer


def support(points, neighbors):
    return frozenset(points) | frozenset(w for v in points for w in neighbors[v])


def connected_word(prefix, points, neighbors):
    assert set(prefix) <= set(points)
    seen = set(prefix)
    out = list(prefix)
    queue = deque(prefix)
    while queue:
        v = queue.popleft()
        new = sorted((neighbors[v] & points) - seen)
        seen.update(new)
        out.extend(new)
        queue.extend(new)
    assert seen == points
    return out


def allowed_count(points, constrained, neighbors, initial=()):
    seq = sorted(points)
    ids = {v:i for i,v in enumerate(seq)}
    masks = [sum(1<<ids[w] for w in neighbors[v] & points) for v in seq]
    optional = [v not in constrained for v in seq]
    seen = sum(1<<ids[v] for v in initial)
    full = (1<<len(seq))-1
    dp = [0]*(full+1)
    dp[seen] = 1
    for mask,value in enumerate(dp):
        if value:
            todo = full ^ mask
            while todo:
                bit = todo & -todo
                todo ^= bit
                i = bit.bit_length()-1
                if optional[i] or masks[i] & mask:
                    dp[mask|bit] += value
    return dp[full]


def completion_h(vertices,seed,neighbors):
    return allowed_count(vertices,vertices,neighbors,seed)


def local_nopeaks(constraints,neighbors):
    points = support(constraints,neighbors)
    return Fraction(allowed_count(points,constraints,neighbors), factorial(len(points)))


def dependency(neighbors):
    return {v:support(neighbors[v] | {v},neighbors)-{v} for v in neighbors}


def normalized_completion(seed,vertices,neighbors,adj2):
    U = vertices-support(seed,neighbors)
    value = Fraction(1)
    for C in components(U,adj2):
        factor = local_nopeaks(C,neighbors)
        assert factor >= Fraction(1,factorial(len(C)+1))
        value *= factor
    return value


def audit_completion():
    total = 0
    for n in (2,3,4):
        vertices,neighbors = square(n)
        adj2 = dependency(neighbors)
        ordered = sorted(vertices)
        if n<=3:
            masks = list(range(1,1<<len(ordered)))
        else:
            rng=Random(20261118)
            masks = [rng.randrange(1,1<<len(ordered)) for _ in range(100)]
        ratios={}
        for mask in masks:
            S=frozenset(v for i,v in enumerate(ordered) if mask&(1<<i))
            expected = Fraction(completion_h(vertices,S,neighbors),factorial(len(vertices)-len(S)))
            actual = normalized_completion(S,vertices,neighbors,adj2)
            assert actual == expected
            ratios[S]=actual
            total+=1
        tested = list(ratios)
        for i in range(len(tested)):
            S,T=tested[i],tested[-i-1]
            U,V=vertices-support(S,neighbors),vertices-support(T,neighbors)
            affected=[C for C in components(U|V,adj2) if C&(U^V)]
            M=sum(map(len,affected))
            assert Fraction(1,factorial(M+1))<=ratios[S]/ratios[T]<=factorial(M+1)
    return total


def hole_pattern(kind):
    _, neighbors=square(12)
    patterns={
        "independent":({(3,3),(8,8)},set(),Fraction(16,25)),
        "axial":({(3,3),(5,3)},{(8,8)},Fraction(29,45)),
        "diagonal":({(3,3),(4,4)},{(8,3),(8,8)},Fraction(13,20)),
        "adjacent":({(3,3),(4,3)},{(8,3),(8,8)},Fraction(3,5)),
    }
    U,inert,factor=patterns[kind]
    holes=support(U,neighbors)|inert
    assert len(holes)==10
    assert local_nopeaks(U,neighbors)==factor
    return U,holes,factor


def audit_tiles():
    results=[]
    for L in (1,2,3,6):
        m,n=12*L,24*L+1
        vertices,neighbors=square(n)
        dep=dependency(neighbors)
        digest={}
        for kind in ("independent","axial","diagonal","adjacent"):
            Ulocal, pattern, factor=hole_pattern(kind)
            W=(m,12*L)
            frames=[(base+12*i,12*j) for base in (0,m+1) for i in range(L) for j in range(2*L)]
            holes={(ox+x,oy+y) for ox,oy in frames for x,y in pattern}
            Utarget={(ox+x,oy+y) for ox,oy in frames for x,y in Ulocal}
            A=frozenset(v for v in vertices-holes if v[0]<m)
            B=frozenset(v for v in vertices-holes if v[0]>m)
            S=A|B|{W}
            assert len(A)==len(B)==m*n-20*L*L
            assert not any(neighbors[v]&B for v in A)
            assert neighbors[W]&A and neighbors[W]&B
            assert vertices-support(S,neighbors)==Utarget
            assert len(vertices-S)==40*L*L+n-1
            p, q=(m-1,12*L-1),(m+1,12*L)
            q_neighbor=(m+1,12*L+1)
            aa=connected_word([p],A,neighbors)
            bb=connected_word([q,q_neighbor],B,neighbors)
            prefix=[p,q,q_neighbor]+aa[1:]+bb[2:]+[W]
            word=connected_word(prefix,vertices,neighbors)
            rank={v:i for i,v in enumerate(word,1)}
            assert {v for v in vertices if all(rank[v]<rank[y] for y in neighbors[v])}=={p,q}
            assert min(rank[y] for y in neighbors[q])==3 and rank[W]==len(S)
            assert rank[p]==1 and rank[q]==2
            for branch in A,B:
                low_degree={v for v in branch if len(neighbors[v]&branch)==2}
                assert all(len(C)<=144 for C in components(low_degree,neighbors))
            Cs=components(Utarget,dep)
            nf=Fraction(1)
            for C in Cs:
                pts=support(C,neighbors)
                assert pts <= holes
                nf*=local_nopeaks(C,neighbors)
            assert nf==factor**(4*L*L)
            digest[kind]=nf
        assert digest["diagonal"]/digest["axial"] == Fraction(117,116)**(4*L*L)
        results.append((n,len(frames),2*len(frames),len(vertices-S)))
    return results


if __name__=="__main__":
    for kind in ("independent","axial","diagonal","adjacent"):
        print(kind,hole_pattern(kind)[2])
    print("exact small completion ratios and factorization cases:",audit_completion())
    print("thick paired states (n, tiles, undominated, suffix):",audit_tiles())
