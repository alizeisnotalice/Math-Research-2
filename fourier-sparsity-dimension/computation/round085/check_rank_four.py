"""python3 -B check_rank_four.py. No third-party dependencies."""
from fractions import Fraction
from itertools import combinations
import json
import platform
from check_components import components, membership, subspaces, vector_rank, wedge


def reduce(x,rows):
    for r in reversed(rows):
        if x & (1 << (r.bit_length()-1)):
            x ^= r
    return x


def echelon(values):
    pivot = {}
    for x in values:
        while x and x.bit_length() in pivot:
            x ^= pivot[x.bit_length()]
        if x:
            pivot[x.bit_length()] = x
    return [pivot[k] for k in sorted(pivot)]


def span(xs):
    out=[0]
    for x in xs:
        out+=[v^x for v in out]
    return out


def pack(columns,n):
    return sum(x<<(n*j) for j,x in enumerate(columns))


def unpack(matrix,n):
    return [(matrix>>(n*j))&((1<<n)-1) for j in range(n)]


def check():
    n=6
    labels = [[1<<(j%n) for j in range(i,i+4)] for i in range(6)]
    labels.extend([[3,6,12,24],[33,18,12,5]])
    N=len(labels)
    m=N.bit_length()-1
    assert all(vector_rank(slot)==4 for slot in labels)
    planes=[span(xs) for xs in labels]
    forms=[wedge(a,b,n)^wedge(c,d,n) for a,b,c,d in labels]
    old=max(Fraction((N+m-q-j+1)*N,N*N-q*(1<<j))
            for q in range(N+1) for j in range(m+1) if q<N or j<m)
    worst=Fraction(0)
    cases=0
    Hcases=0
    for rows in subspaces(n):
        proj=[[reduce(column,rows) for column in unpack(mat,n)] for mat in forms]
        ps=[vector_rank(cols) for cols in proj]
        k=vector_rank(pack(cols,n) for cols in proj)
        low=[cols for cols,p in zip(proj,ps) if p==1]
        k1=vector_rank(pack(cols,n) for cols in low)
        info=echelon(x for cols in low for x in cols)
        cut=[pack([reduce(x,info) for x in cols],n) for cols in proj]
        pieces=components([x for x in cut if x])
        assert k >= k1+sum(len(part) for part in pieces)
        for part in pieces:
            assert vector_rank(x for w in part for x in unpack(w,n))<=3*len(part)+1
        t=n-len(rows)
        assert t<=k1+sum(3*len(part)+1 for part in pieces)
        if t:
            assert k+1+2*sum(len(part)+1 for part in pieces)>=t+1
        hits=sorted((sum(membership(x,rows) for x in S) for S in planes),reverse=True)
        count=[0]
        for h in hits:
            count.append(count[-1]+h)
        for j in range(m+1):
            for q,inside in enumerate(count):
                codim=N+m+n-q-j-len(rows)
                if not codim:
                    continue
                outside=16*N*N-(1<<j)*inside
                num=(codim+1)*4*N
                assert num*old.denominator<=old.numerator*outside
                if num*worst.denominator>worst.numerator*outside:
                    worst=Fraction(num,outside)
                cases+=1
        Hcases+=1
    return dict(python=platform.python_version(), seed=None,tolerance=0,
                scope='actual pure rank-four address-selector with 8 independent guards and 6 phase coordinates',
                old_L=str(old),new_L=str(worst),address_kernel_checks=Hcases,
                exact_reduced_cut_checks=cases)


if __name__=='__main__':
    print(json.dumps(check(),indent=2))
