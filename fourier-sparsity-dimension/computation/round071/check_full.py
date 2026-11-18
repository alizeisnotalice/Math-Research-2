"""Check all coordinate splits of all four-bit truth tables exactly."""
from collections import Counter
from functools import cache
import json
import platform


@cache
def spectral(table, size):
    row = [1-2*((table>>i)&1) for i in range(size)]
    width = 1
    while width < size:
        for i in range(0, size, 2*width):
            for j in range(i,i+width):
                a,b=row[j],row[j+width]
                row[j],row[j+width]=a+b,a-b
        width *= 2
    return tuple(i for i,x in enumerate(row) if x)


def rank(xs):
    pivots={}
    for x in xs:
        while x and x.bit_length() in pivots:
            x ^= pivots[x.bit_length()]
        if x:
            pivots[x.bit_length()]=x
    return len(pivots)


@cache
def filled(labels, size):
    groups={}
    for b in labels:
        groups.setdefault(b%size,[]).append(b//size)
    return sum(1 << rank(v^group[0] for v in group) for group in groups.values())


def main():
    counts=Counter()
    for t in range(1<<16):
        labels=spectral(t,16)
        d=rank(labels)
        for k in range(5):
            size=1<<k
            L=max(len(spectral((t>>i)&((1<<size)-1),size))
                  for i in range(0,16,size))
            M=filled(labels,size)
            assert len(labels)<=M and d*d<=7*L*M
            counts['splits']+=1
            counts['with_holes']+=int(M>len(labels))
            counts['character_only']+=int(L==1)
    return dict(python=platform.python_version(), n=4,tables=65536,**counts)


if __name__=='__main__':
    print(json.dumps(main(),indent=2))
