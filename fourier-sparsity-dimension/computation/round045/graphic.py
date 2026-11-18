"""Enumerate all bonds of necklaces with three through eight triangles, modulo cut complement."""
import json
import platform

def check(m):
    edges=[]
    for i in range(m):
        j=(i+1)%m
        edges.extend(((i,j),(i,m+i),(j,m+i)))
    n=2*m
    mask=(1<<n)-1
    adj=[0]*n
    for i,j in edges:
        adj[i]|=1<<j
        adj[j]|=1<<i
    def connected(side):
        seen=front=side & -side
        while front:
            expand=0
            while front:
                bit=front & -front
                expand |= adj[bit.bit_length()-1]
                front^=bit
            front=expand & side & ~seen
            seen|=front
        return seen==side
    counts={}
    for left in range(1,mask,2):
        if connected(left) and connected(mask^left):
            cut=sum(((left>>i)^(left>>j))&1 for i,j in edges)
            counts[cut]=counts.get(cut,0)+1
    assert counts=={2:m,4:2*m*(m-1)}
    return dict(triangles=m,rank=n-1,nullity=len(edges)-n+1,bonds=counts)

print(json.dumps(dict(python=platform.python_version(),graphs=[check(m) for m in range(3,9)]),indent=2))
