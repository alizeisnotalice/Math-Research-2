"""Exact three-valued support checks and examples with a minimal odd core."""
from collections import Counter
from itertools import product
import json
import platform


def transform(row):
    vals=list(row)
    s=1
    while s<len(vals):
        for j in range(0,len(vals),s*2):
            for i in range(j,j+s):
                a,b=vals[i],vals[i+s]
                vals[i],vals[i+s]=a+b,a-b
        s*=2
    return vals


def example(a,normals,complement):
    physical=1<<a
    h=[int(((u%(1<<normals))==0) != complement) for u in range(physical)]
    g=[1-v for v in h]
    g[g.index(1)]=-1
    truth=[g[u]+(1-2*x)*h[u] for x in (0,1) for u in range(physical)]
    assert set(truth)=={-1,1}
    w=transform(truth)
    assert all(v%4==0 for v in w)
    odd={i for i,v in enumerate(w) if v//4%2}
    tail=sum(bool(v) for v in w[physical:])
    assert odd==set(range(physical)) and tail==1<<normals
    return dict(a=a,normal_rank=normals,complement=complement,core=len(odd),tail=tail)


counts=Counter()
for row in product((-1,0,1),repeat=8):
    s=sum(bool(v) for v in transform(row))
    counts[s]+=1
    if s==2:
        assert sum(bool(v) for v in row)==4
    if s==4:
        assert sum(bool(v) for v in row) in (2,6,8)
assert counts[3]==0
print(json.dumps(dict(python=platform.python_version(),ternary_3var=dict(sorted(counts.items())),
    examples=[example(a,b,c) for a in (4,5,6) for b,c in [(1,False),(2,False),(2,True)]]),indent=2))
