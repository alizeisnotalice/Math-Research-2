"""Verify the common-kernel graph example with a full integer Walsh transform."""
from itertools import product
import json
import platform


def walsh(values):
    size = 1
    while size < len(values):
        for i in range(0, len(values), 2*size):
            for j in range(i,i+size):
                a,b=values[j],values[j+size]
                values[j],values[j+size]=a+b,a-b
        size *= 2
    return values


def rank(vectors):
    basis = {}
    for v in vectors:
        while v:
            j = v.bit_length()
            if j not in basis:
                basis[j] = v
                break
            v ^= basis[j]
    return len(basis)


def model(b,p):
    sizes=[2**b,2**b,2**p,2,2,2**b,2**p]
    values=[]
    for u,v,h,y1,y2,z,t in product(*map(range,sizes)):
        inside=(t==h and (u&z).bit_count()%2==y1 and (v&z).bit_count()%2==y2)
        values.append(1-2*inside)
    width=2**(2+b+p)
    bands={}
    for i in range(0,len(values),width):
        w=walsh(values[i:i+width])
        assert sum(map(abs,w))<3*width
        for j,c in enumerate(w):
            if c:
                bands.setdefault(j,set()).add(i//width)
    domains=set(map(frozenset,bands.values()))
    hits={sum(j in a for a in domains) for j in range(len(values)//width)}
    s=[j for j,w in enumerate(walsh(values)) if w]
    assert hits=={4} and len(domains)==1+3*2**b
    assert (len(s),rank(s))==(2**p*(1+3*4**b),3*b+2+p)
    return dict(b=b,p=p,domains=len(domains),active_domains_per_query=4,s=len(s),d=rank(s))


print(json.dumps(dict(python=platform.python_version(), models=[
    model(b,p) for b,p in [(1,1),(1,2),(2,1),(2,2),(3,2)]]),indent=2))
