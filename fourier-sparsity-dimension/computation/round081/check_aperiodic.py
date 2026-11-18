"""Exact connectivity checks and Walsh replay of the period-free selector."""
from collections import Counter
from itertools import combinations
import json
import platform


def rank(xs):
    basis={}
    for x in xs:
        while x and x.bit_length() in basis:
            x^=basis[x.bit_length()]
        if x:
            basis[x.bit_length()]=x
    return len(basis)


def walsh(f):
    w=list(f)
    step=1
    while step<len(w):
        for i in range(0,len(w),2*step):
            for j in range(i,i+step):
                a,b=w[j],w[j+step]
                w[j],w[j+step]=a+b,a-b
        step*=2
    return w


def spaces(n):
    last={frozenset([0])}
    full=list(last)
    for _ in range(n):
        last={frozenset(H|{p^a for a in H})
              for H in last for p in range(1 << n) if p not in H}
        full.extend(last)
    return full


def rectangle_connected(X):
    parts={x:{x} for x in X}
    for t in combinations(X,4):
        if t[0]^t[1]^t[2]^t[3]:
            continue
        grp=set().union(*(parts[x] for x in t))
        for x in grp:
            parts[x]=grp
    return len(parts[next(iter(X))])==len(X)


def small_all():
    kernels=spaces(4)
    count=0
    for mask in range(16):
        if mask.bit_count()<2:
            continue
        S={x | b<<2 for b in range(4) if mask>>b&1 for x in range(4)}
        for H in kernels:
            X=S-H
            if X:
                assert rectangle_connected(X)
                count+=1
    return dict(base_sets=11,kernels=len(kernels),nonempty_exteriors=count)


def all_colorings():
    n=4
    N=1<<n
    cases={}
    for cut in (0,1,2,4,8):
        X=list(range(cut,N))
        size=len(X)
        R=rank([N|x for x in X])
        checks=0
        best=1000
        # Complementary colorings are identical for this symmetric question.
        for z in range(1<<(size-1)):
            Y=[N|x for i,x in enumerate(X) if z>>i&1]
            Z=[N|x for i,x in enumerate(X) if not(z>>i&1)]
            y,zr=rank(Y),rank(Z)
            if max(y,zr)<R:
                lam=y+zr-R
                assert lam>=2
                best=min(best,lam)
                checks+=1
        cases[cut]=dict(essential=checks,min_intersection=best)
    return cases


def gadget(m):
    M=1<<m
    A=M*M
    f=[1-2*(((x&(M-1))&(x>>m)).bit_count()%2) for x in range(A)]
    g=[1-2*(x < M//2) for x in range(A)]
    table=f+g
    spec=walsh(table)
    T={x for x,v in enumerate(spec) if v}
    assert 0 in T and 2*A-len(T)==2*M-1
    missed=set(range(2*A))-T
    assert not any(all(x^p in missed for x in missed) for p in range(1,2*A))
    D=[int(t<0) for t in table]
    assert {i for i,v in enumerate(walsh(D)) if v}==T
    assert {i for i,v in enumerate(walsh([1-v for v in D])) if v}==T
    return D,T


def real_selector(m):
    D,T=gadget(m)
    K=len(D)
    # On three base data bits, a single-point flip is odd and full spectrum.
    F=[-1]+[1]*7
    src={x for x,c in enumerate(walsh(F)) if c}
    parent=[(F[x&7] if D[x>>4] else 1-2*((x>>3)&1)) for x in range(16*K)]
    W=walsh(parent)
    S={x for x,c in enumerate(W) if c}
    assert S=={x|(t<<4) for x in src|{8} for t in T}
    assert rank(S)==rank(src)+1+2*m+1
    assert len(S)==(len(src)+1)*(K-(1<<(m+1))+1)
    assert all(parent[x|(y<<4)]==F[x] for y in range(K) if D[y] for x in range(8))
    return dict(k=2*m+1,source_s=len(src),outer_factor_s=len(T),
                selector_s=len(S),selector_d=rank(S),full_truths=len(parent))


if __name__=='__main__':
    print(json.dumps(dict(python=platform.python_version(),tolerance=0,
                          four_times_kernel_test=small_all(),
                          all_local_unordered_two_colorings=all_colorings(),
                          period_free_selectors=[real_selector(m) for m in (3,4,5)]),indent=2))
