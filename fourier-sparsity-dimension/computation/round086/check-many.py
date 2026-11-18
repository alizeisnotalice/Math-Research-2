from fractions import Fraction
from itertools import combinations
from random import Random
import json
import platform


def reduce(v, piv):
    while v:
        j=v.bit_length()
        if j not in piv:
            return v
        v ^= piv[j]
    return v


def basis(rows):
    piv={}
    for v in rows:
        v=reduce(v,piv)
        if v:
            piv[v.bit_length()]=v
    return piv


def check():
    rng=Random(86052)
    edges=[(i,(i+1)%5) for i in range(5)]+[(5+i,5+(i+2)%5) for i in range(5)]+[(i,5+i) for i in range(5)]
    edge_masks=[(1<<i)|(1<<j) for i,j in edges]
    assert min(sum((x&m).bit_count()%2 for m in edge_masks) for x in range(1,1023))==3
    counts={}
    for r in (10,11):
        hist=[0]*15
        cnt=0
        for ix in combinations(range(15),r):
            cols=[edge_masks[i] for i in ix]
            if all(len(basis(cols[:j]+cols[j+1:]))==9 for j in range(r)):
                cnt+=1
                for i in ix:
                    hist[i]+=1
        assert len(set(hist))==1
        counts[r]=(cnt,hist[0])
    assert counts[10][0]==0 and Fraction(counts[11][1],counts[11][0])==Fraction(11,15)
    params=[]
    tested=0
    for m in (1,2,3,7,17):
        a=(15*m).bit_length()
        p=(1<<a)-15*m
        for shared in (False,True):
            n=2
            def R(n):
                return 9*m+a+(n+p if shared else p*(n+1))
            while (1<<(n//2)) < 3*(R(n)+1):
                n+=2
            assert n<=4*a+8
            t=a+n
            inside=[0,1]+[v for j in range(1,t) for v in (1<<j,1+(1<<j))]
            lifted=[(1<<t)+v for v in inside]
            assert len(inside)==2*t and all(len(basis(lifted[:j]+lifted[j+1:]))==t+1 for j in range(2*t))
            assert 2*t*15<=11*(1<<(n//2))
            params.append(dict(m=m,a=a,p=p,shared=shared,minimum_even_n=n,R=R(n),fill_pick=2*t))
        # Random kernels of maps out of the direct sum, including mixed-block maps.
        # Coordinate j is represented by e_j+e_9 in each 9-dimensional edge space.
        for cod in range(1,min(9*m,15)+1):
            for _ in range(14):
                images=[[rng.randrange(1<<cod) for j in range(9)]+[0] for i in range(m)]
                global_rank=len(basis(x for row in images for x in row))
                if not global_rank:
                    continue
                mass=sum(row[i]!=row[j] for row in images for i,j in edges)
                assert 3*(global_rank+1)<=2*mass
                tested+=1
    return dict(python=platform.python_version(),seed=86052,graphic=counts,quotient_cuts=tested,parameters=params)


if __name__=='__main__':
    print(json.dumps(check(),indent=2))
