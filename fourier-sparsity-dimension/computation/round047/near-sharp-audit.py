"""Exact finite verification of the replacement, rank and obstruction."""
import json
import platform
import random

def fwht(values):
    out=values.copy()
    step=1
    while step<len(out):
        for start in range(0,len(out),2*step):
            for j in range(start,start+step):
                out[j],out[j+step]=out[j]+out[j+step],out[j]-out[j+step]
        step*=2
    return out

def f2rank(X):
    basis={}
    for a in X:
        while a and a.bit_length() in basis:
            a^=basis[a.bit_length()]
        if a: basis[a.bit_length()]=a
    return len(basis)

def linear_spaces(n,k):
    spaces={frozenset({0})}
    for dim in range(k):
        spaces={H|{x^v for x in H} for H in spaces for v in range(1,1<<n) if v not in H}
    return spaces

n=6
rng=random.Random(4701)
g=[1-2*rng.randrange(2) for x in range(1<<n)]
max_monochromatic=0
flat_census={}
for k in range(5):
    c=0
    for H in linear_spaces(n,k):
        pool=set(range(1<<n))
        while pool:
            x=next(iter(pool))
            A={x^s for s in H}
            pool-=A
            if len({g[v] for v in A})==1: c+=1
    flat_census[k]=c
    if c:max_monochromatic=k
assert max_monochromatic<4
G=[i for i,c in enumerate(fwht(g)) if c]
# Height 3 ; delete whole left depth-one subtree. Outside internal
# coordinates are root, right child, two right grandchildren (=bits 0..3).
vals=[]
for x in range(1<<(n+4)):
    z=x&15
    y=x>>4
    vals.append(g[y] if not z&1 else 1-2*((z>>(2+((z>>1)&1)))&1))
S=[i for i,c in enumerate(fwht(vals)) if c]
assert len(S)==8+2*len(G)
assert f2rank(S)==4+f2rank(G)
print(json.dumps(dict(python=platform.python_version(),seed=4701,n=n,height=3,cut_depth=1,
        walsh_sparsity_g=len(G),walsh_rank_g=f2rank(G),walsh_sparsity_F=len(S),walsh_rank_F=f2rank(S),
        monochromatic_flat_counts=flat_census,maximum_monochromatic_dimension=max_monochromatic,
        lower_bound_constant_leaves=2**(n-max_monochromatic)),indent=2))
