from collections import defaultdict, deque
from fractions import Fraction
from itertools import permutations
from math import comb


def masks_grid(n):
    coords = list((x,y) for y in range(n) for x in range(n))
    index = {p:i for i,p in enumerate(coords)}
    masks = []
    for x,y in coords:
        masks.append(sum(1<<index[q] for q in ((x-1,y),(x+1,y),(x,y-1),(x,y+1)) if q in index))
    return coords,masks


def fibre_counts(n):
    coords,masks = masks_grid(n)
    groups=defaultdict(dict)
    for word in permutations(range(n*n)):
        seen=0
        roots=[]
        merge=None
        comps=[]
        for i,v in enumerate(word):
            mask=masks[v]
            if not (seen&mask):
                roots.append(v)
                if len(roots)>2:
                    break
            seen |= 1<<v
            if merge is None:
                touches=[c for c in comps if c&mask]
                if len(touches)>=2:
                    merge=(i+1,v,tuple(sorted(touches)))
                else:
                    if touches:
                        comps.remove(touches[0])
                    comps.append((1<<v)|(touches[0] if touches else 0))
        if len(roots)!=2:
            continue
        p,q=roots
        rank=word.index(q)+1
        key=q,tuple(v for v in word if v!=q)
        distance=sum(abs(x-y) for x,y in zip(coords[p],coords[q]))
        neighbor=next(i for i,v in enumerate(word,1) if masks[q]&(1<<v))
        groups[key][rank]=(distance,neighbor,merge)
    count=0
    for fiber in groups.values():
        values=set(fiber.values())
        assert len(values)==1
        d,T,_=values.pop()
        assert set(fiber)==set(range(2,T))
        L=T-2
        moment=sum((rank-1)**(d-1) for rank in fiber)
        assert d*moment>=L**d
        for s in range(2,len(coords)+1):
            early=sum(rank<=s for rank in fiber)
            assert d*(s-1)*moment>=L**d*early
        count+=len(fiber)
    return len(groups),count


def binomial_checks():
    count=0
    for m in range(2,81):
        for b in range(1,m):
            pmf={
                j:Fraction(comb(m-1-j,b-1),comb(m-1,b))
                for j in range(1,m-b+1)
            }
            assert sum(pmf.values())==1
            ell=(m-1)//(2*b)+1
            assert sum(prob for j,prob in pmf.items() if j>=ell)>=Fraction(1,2)
            for d in (2,3,5):
                moment=sum(j**(d-1)*prob for j,prob in pmf.items())
                for s in (2,3,6,max(2,m//2),m):
                    early=sum(prob for j,prob in pmf.items() if j<=s-1)
                    assert moment>=early*Fraction(m,b)**d/(4**d*(s-1))
                    count+=1
    return count


def comb_tree(L):
    arm_rows=range(2*L+2,3*L+1,2)
    return (
        {(0,y) for y in range(-3*L,3*L+1)}
        | {(-x,0) for x in range(1,3*L+1)}
        | {(-x,sgn*y) for y in arm_rows for x in range(1,L+1) for sgn in (-1,1)}
        | {(-x,y) for x in arm_rows for y in range(1,L+1)}
    )


def local_neighbors(vertices,v):
    x,y=v
    return vertices & {(x-1,y),(x+1,y),(x,y-1),(x,y+1)}


def rooted_tree_ratio(vertices,root):
    neighbors={v:local_neighbors(vertices,v) for v in vertices}
    parent={root:None}
    depth={root:0}
    order=[root]
    for v in order:
        for w in neighbors[v]:
            if w not in parent:
                parent[w]=v
                depth[w]=depth[v]+1
                order.append(w)
    assert len(order)==len(vertices)
    assert sum(map(len,neighbors.values()))==2*(len(vertices)-1)
    size={v:1 for v in vertices}
    for v in order[:0:-1]:
        size[parent[v]]+=size[v]
    ratio={root:Fraction(1)}
    for v in order[1:]:
        ratio[v]=ratio[parent[v]]*Fraction(size[v],len(vertices)-size[v])
        assert ratio[v]<=Fraction(1,2)**depth[v]
    normalizer=sum(ratio.values())
    assert normalizer<=9
    return normalizer


def audit_combs():
    entries=[]
    for d in (3,4,8):
        for L in (8,16,32,64):
            core=comb_tree(L)
            a=core|{(x,0) for x in range(1,d-1)}
            b={(d-x,y) for x,y in core}
            bridge=(d-1,0)
            assert not (a&b)
            assert all(not local_neighbors(b,v) for v in a)
            assert len(local_neighbors(a,bridge))==len(local_neighbors(b,bridge))==1
            t=3*L+L*(L//2)
            assert len(a)==1+3*t+d-2 and len(b)==1+3*t
            normalizer_a=rooted_tree_ratio(a,(0,0))
            normalizer_b=rooted_tree_ratio(b,(d,0))
            a_n,b_n=len(a),len(b)
            m=a_n+b_n
            theta=Fraction(2*a_n*b_n,m*(m-1))
            assert theta/(normalizer_a*normalizer_b)>=theta/81
            entries.append((6*L+d+1,d,m,float(theta/(normalizer_a*normalizer_b))))
    return entries


if __name__=="__main__":
    for n in (2,3):
        print(f"n={n}: exact parent fibres, original K2 words = {fibre_counts(n)}")
    print("without-replacement moment inequality checks:",binomial_checks())
    for entry in audit_combs():
        print("n, centre distance, premerge size, exact-fraction near-top2 probability (display):", entry)
