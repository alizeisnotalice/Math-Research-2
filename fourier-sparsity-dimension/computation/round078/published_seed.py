"""Example 2, P_0 of Wang–Zhang–Tang–Zhou v2; integer transform and density."""
from deletion import rank,spaces
from collections import Counter
from fractions import Fraction
import json

def fwht(v):
    v=list(v)
    k=1
    while k<len(v):
        for i in range(0,len(v),2*k):
            for j in range(i,i+k):
                v[j],v[j+k]=v[j]+v[j+k],v[j]-v[j+k]
        k*=2
    return v

def main():
    p=4
    def label(u):
        return u + ((u&3==3)<<p) + ((u==15)<<(p+1))
    labels={label(u) for u in range(16)}
    f=[(-1)**((label(i&15)&(i>>4)).bit_count()%2) for i in range(1024)]
    w=fwht(f)
    s=[i for i,t in enumerate(w) if t]
    assert all(w[(label(u)<<4)+a]==64*(-1)**((a&u).bit_count()%2) for u in range(16) for a in range(16))
    assert len(s)==256 and rank(s)==10
    assert all(sum(w[a]*w[a^x] for a in s)==(1024**2 if x==0 else 0) for x in range(1024))
    anf=[int(v<0) for v in f]
    for bit in range(10):
        for m in range(1024):
            if m>>bit&1: anf[m]^=anf[m^(1<<bit)]
    density,u,j=max((Fraction((1<<j)*len(labels&space)-1,j+rank(space)),tuple(sorted(space)),j)
                    for space in spaces(6) for j in range(5) if j+rank(space))
    print(json.dumps(dict(seed="Wang et al. v2 Example 2 P0",s=len(s),d=rank(s),
        degree=max(m.bit_count() for m,v in enumerate(anf) if v), spectrum_hull_rank=rank(i^s[0] for i in s),
        signed_integer_walsh={str(k):v for k,v in Counter(w).items()}, coefficients=w,
        delta_star="15/32",density=str(density),c=-(-density.numerator//density.denominator),
        core_dimension=j,dense_label_space=u,nonzero_convolutions_checked=1023,
        K=4, relation_dimension=0, intersection_dimension=0,lift_d=20,lift_c=13),indent=2))
if __name__=='__main__': main()
