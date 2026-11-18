"""Integer verification of a quotient band separated in the augmented matroid."""
from collections import Counter
from functools import cache
import json
import platform


def walsh(values):
    w=list(values)
    k=1
    while k<len(w):
        for i in range(0,len(w),2*k):
            for j in range(i,i+k):
                lo,hi=w[j],w[j+k]
                w[j],w[j+k]=lo+hi,lo-hi
        k*=2
    return w


def rank(vec):
    piv={}
    for v in vec:
        while v and v.bit_length() in piv:
            v ^= piv[v.bit_length()]
        if v:
            piv[v.bit_length()]=v
    return len(piv)


@cache
def setup(sup,N,k):
    X=[a for a in sup if a >> k]
    aug=lambda S:rank(N|a for a in S)
    rx=aug(X)
    r=rank(a >> k for a in X)
    ans=[]
    if r>=2:
        for beta in set(a >> k for a in X):
            B=tuple(a for a in X if a>>k==beta)
            C=tuple(a for a in X if a>>k!=beta)
            if aug(B)+aug(C)==rx:
                ans.append((beta,B,C))
    if r>=3:
        assert len(ans)<=1
    return r,ans


def check(f,k,validate_all=True,target=None):
    N=len(f)
    w=walsh(f)
    S=tuple(a for a,c in enumerate(w) if c)
    r,ans=setup(S,N,k)
    out=[]
    for beta,B,C in ans:
        if target is not None and beta!=target:
            continue
        v=walsh([2*c if a in B else 0 for a,c in enumerate(w)])
        p=walsh([c if a in C else 0 for a,c in enumerate(w)])
        assert all(abs(i)==N for i in v) and all(i in (-N,0,N) for i in p)
        q=[1-2*(t//N)**2 for t in p]
        wq=walsh(q)
        T=[i for i,c in enumerate(wq) if c]
        assert all(i >> k == beta for i in T)
        assert set(S)==set(B)|set(C)|{b^t for b in B for t in T}
        assert len(S)==len(C)+len(B)*(1+len(T))
        assert all(2*N*z==2*y+x*(1+a) for z,y,x,a in zip(f,p,v,q))
        if r>=3:
            assert all(2*abs(w[a])<N for a in C)
            if len(B)>1:
                assert all(2*abs(c)!=N for c in w)
        b_list=B if validate_all else B[:1]
        for b in b_list:
            f0=[y//N + (1-2*((b&i).bit_count()%2))*(1+a)//2
                for i,(a,y) in enumerate(zip(q,p))]
            assert set(f0)<={1,-1}
            w0=walsh(f0)
            S0=[i for i,c in enumerate(w0) if c]
            assert set(S0)==set(C)|{b}|{b^t for t in T}
            assert 2*abs(w0[b])==N
            aff=lambda s:rank(N|a for a in s)-1
            assert aff(S)-aff(S0)==aff(B)
            assert rank(S)-rank(S0) in (aff(B),aff(B)-1)
        out.append((r,len(B),len(T),len(S),rank(S)))
    return out


def exhaust_four():
    stats=Counter()
    for table in range(1 << 16):
        f=[1-2*((table >> i)&1) for i in range(16)]
        for k in (1,2):
            for r,b,t,s,d in check(f,k):
                stats[(k,r,b,t)]+=1
    return {','.join(map(str,key)):val for key,val in sorted(stats.items())}


def constructed():
    # Kernel w-input, E-input, and the private selector u. The first high
    # frequency is t and the remaining two high frequencies are z1,z2.
    seeds=[
        (1,lambda x:1-2*x),
        (2,lambda x:1-2*(x==0)),
        (3,lambda x:1-2*(x==0)),
        (4,lambda x:1-2*(((x>>3)&1)*(((x>>2)&1)^((x&1)*((x>>1)&1))))),
    ]
    tests=[]
    for nw,W in seeds:
        for ne in (0,2):
            k=nw+ne+1
            f=[]
            for x in range(1 << (k+3)):
                w=W(x&((1<<nw)-1))
                e=1-2*int(ne!=0 and (x>>nw)&((1<<ne)-1)==0)
                u=(x>>(nw+ne))&1
                t=1-2*((x>>k)&1)
                g=1-2*((x>>(k+1+u))&1)
                f.append(((e+t)*w+(1-t*e)*g)//2)
            row,=check(f,k,target=1)
            tests.append(row)
    return tests


if __name__=='__main__':
    print(json.dumps(dict(python=platform.python_version(),tolerance=0,
                          four_bit_tables=65536,fixed_kernel_dimensions=[1,2],
                          accepted_by_kernel_quotient_band_T=exhaust_four(),
                          examples_columns=['r','band','T','s','d'],
                          examples=constructed()),indent=2))
