"""Exact rational and permutation checks for the exclusive witness bound."""
from collections import Counter
from fractions import Fraction as F
from itertools import permutations, product
from math import comb


def fiber():
    points=[F(i,18) for i in range(18)]
    for a,b,c in product(points,repeat=3):
        u,v=1-a,1-b
        if b<=a:
            before=(1-c)*u
            sufficient=(1-c)*u*u/2
        else:
            L=1-max(b,c)
            J=max(b-c,0)*v
            before=L*u+J
            sufficient=L*(u-v+v*v/2)+J*v/2
        assert 2*sufficient>=u*before
    print('rational_fibers', len(points)**3)


def full_words():
    # d,t opposite u,v,w; z shares the constrained-side parity.
    d,t,z,u,v,w,h=range(7)
    edges={(d,z),(d,u),(d,v),(t,v),(t,w),(h,z),(h,u)}
    adj={a:{y if x==a else x for x,y in edges if a in (x,y)} for a in range(7)}
    numerator=denominator=rank_sum=0
    for word in permutations(range(7)):
        ranks={p:i for i,p in enumerate(word,1)}
        forbidden=lambda a:all(ranks[a]<ranks[x] for x in adj[a])
        if any(forbidden(p) for p in (u,v,w,h)):
            continue
        denominator+=1
        numerator+=not forbidden(d)
        witness=ranks[u] if ranks[u]<ranks[h] else 8
        rank_sum+=witness
    assert numerator*16>=rank_sum>0
    print('full_words',denominator,numerator,'exclusive_rank_sum',rank_sum)


def grid_family():
    A={(j,1) for j in range(2,5)}|{(2,2),(4,2)}|{(5,y) for y in range(2,6)}
    w=(4,3)
    B={(3,3)}
    d,z,t=(3,4),(4,4),(4,5)
    K=A|{w}
    for n in range(6,41):
        points=set(product(range(n),repeat=2))
        def boundary(S):
            return {v for x,y in S for v in ((x-1,y),(x+1,y),(x,y-1),(x,y+1)) if v in points-S}
        def connected(S):
            reached={min(S)}
            while new:=boundary(reached)&(S-reached):
                reached.update(new)
            return reached==S
        assert connected(A) and boundary(A).isdisjoint(B)
        assert w in boundary(B)&boundary(A)
        extra=boundary(B)-K-boundary(K)
        safe=boundary(K)-B
        U=points-K-B-safe-extra
        assert extra=={d} and {z,t}<=safe and z in boundary({d})
        old_neighbors=boundary({d})&U
        assert old_neighbors=={(2,4),(3,5)}
        assert old_neighbors&boundary({t})=={(3,5)}
    print('grid_partial_backup',35)


def truncated_group_moment():
    tested=0
    for k in range(1,5):
        for g in range(4):
            for groups in product((2,3,4),repeat=g):
                coeff=Counter({0:1})
                for r in groups:
                    new=coeff.copy()
                    for degree, value in coeff.items():
                        new[degree+r]-=value
                    coeff=new
                for a in (F(0), F(1,4), F(1,2), F(3,4), F(7,8)):
                    moment=F(0)
                    for degree, value in coeff.items():
                        for j in range(degree+1):
                            moment += F(k*value*comb(degree,j),k+j)*a**(degree-j)*(1-a)**(k+j)
                    assert 0<(1-a)**(k+g)/comb(k+g,k)<=moment<=(1-a)**k
                    tested+=1
    print('partial_group_moments',tested)


if __name__=='__main__':
    fiber()
    full_words()
    grid_family()
    truncated_group_moment()
