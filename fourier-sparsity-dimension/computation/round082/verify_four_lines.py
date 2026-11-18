"""Exact all-query optimization on recursive four-line selectors."""
from collections import Counter
from fractions import Fraction
from itertools import combinations
import json
import platform

LINES = ((0,1),(2,4),(3,7),(5,6))


def rank(values):
    b={}
    for x in values:
        while x and x.bit_length() in b:
            x^=b[x.bit_length()]
        if x:
            b[x.bit_length()]=x
    return len(b)


def family(depth, sizes=None):
    fresh=0
    def alloc(k):
        nonlocal fresh
        loc=fresh
        fresh+=k
        return loc
    def leaf():
        return {1<<alloc(1):1},1
    def join(blocks, cells, loc, bits):
        den=(1<<bits)//len(cells[0])*max(m for _,m in blocks)
        S={}
        for (sp,m),cell in zip(blocks,cells):
            for a in range(1<<bits):
                c=sum(1-2*((a&x).bit_count()%2) for x in cell)
                if c:
                    for f,x in sp.items():
                        S[f|(a<<loc)]=x*c*den//((1<<bits)*m)
        return S,den
    def binary(n):
        if n==1:
            return leaf()
        loc=alloc(1)
        return join([binary(n//2),binary(n//2)],((0,),(1,)),loc,1)
    def macro(k):
        if not k: return leaf()
        loc=alloc(3)
        return join([macro(k-1) for _ in LINES],LINES,loc,3)
    if sizes is None: S,den=macro(depth)
    else:
        loc=alloc(3)
        S,den=join([binary(n) for n in sizes],LINES,loc,3)
    assert sum(c*c for c in S.values())==den*den
    return S,den,fresh


def cut(S,q):
    top=1<<(q.bit_length()-1)
    a,b={},{}
    for x,v in S.items():
        if x&top:
            z=x^q
            if z in S:continue
            xv,yv=0,v
        else:
            z=x
            xv,yv=v,S.get(x^q,0)
        if xv+yv:a[z]=xv+yv
        if xv-yv:b[z]=xv-yv
    return a,b


def audit(depth,sizes=None,all_directions=True):
    S,den,n=family(depth,sizes)
    tag=1<<n
    rng=lambda pts:rank(tag|p for p in pts)
    R=rng(S);s=len(S);D=R+1
    if all_directions:
        directions={a^b for a,b in combinations(S,2)}
    else:
        directions={q for q in range(1,8)}
    no_collision=(1<<(R-1))-1-len(directions) if all_directions else None
    status=Counter()
    pbest=Fraction((D-1)**2,s) if no_collision else Fraction(0)
    best_n=0
    topdata={}
    # The cut modulo q has two child spans whose union has rank R-1.
    for q in directions:
        A,B=cut(S,q)
        r0,r1=rng(A),rng(B)
        s0,s1=len(A),len(B)
        u=r0+r1+2;l=s0+s1;j=r0+r1-(R-1)
        common=A.keys()&B.keys()
        k=len(common);e=sum(abs(A[p])==abs(B[p]) for p in common)
        assert 2*l-2*k-e==s
        passed=2*l*(2*u*j-j*j)>=(2*k+e)*u*u
        improved=max((r0+1)**2*s*s1,(r1+1)**2*s*s0)>=D*D*s0*s1
        status[str((passed,improved))]+=1
        p=max(Fraction((r0+1)**2,s0), Fraction((r1+1)**2,s1))
        if p>pbest: pbest,best_n=p,1
        elif p==pbest: best_n+=1
        if q<8:
            topdata[q]=[s0,r0+1,s1,r1+1,k,j,e,passed,improved]
    if sizes is None:
        N=4**depth
        assert pbest==Fraction((2*N-1)**2,N*N-3*N//4)
        assert status=={'(False, False)':len(directions)}
    pair_cells=Counter()
    if sizes is None:
        for q,t in combinations(range(1,8),2):
            rem=t^q if t & (1<<(q.bit_length()-1)) else t
            for branch in cut(S,q):
                for leaf in cut(branch,rem):
                    ds=rng(leaf)+1
                    assert ds*ds==4*len(leaf)
                    pair_cells[str((len(leaf),ds))]+=1
    if sizes is not None and sizes[1:]==(4*sizes[0],)*3:
        # Root queries 2 and 4 leave bit 1 free, the direction of LINES[0].
        for first in cut(S,2):
            for last in cut(first,4):
                dd=rng(last)+1
                assert dd**2==4*len(last)
                pair_cells[str((len(last),dd))]+=1
        if all_directions:
            assert status.get('(True, True)',0)==0
            assert status['(False, True)']==3
    return dict(s=s,D=D,denominator=den,affine_directions=(1<<(R-1))-1,
                noncolliding=no_collision,examined_counts=status,
                best_child_P_squared=str(pbest),best_collision_directions=best_n,
                root_results=topdata,two_query_root_cells=pair_cells)


if __name__=='__main__':
    print(json.dumps(dict(python=platform.python_version(),tolerance=0,
        symmetric=[audit(k) for k in (1,2)],
        asymmetric=[audit(0,sizes,False) for sizes in ((1,4,1,4),(2,8,2,8))],
        asymmetric_global=[audit(0,(1,4,4,4)),audit(0,(2,8,8,8),False)]),
        indent=2))
