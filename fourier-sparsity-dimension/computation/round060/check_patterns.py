"""Integer verification of support-pattern atoms, query rank and sparsity."""
import json
import platform
from collections import Counter
from itertools import combinations


def walsh(values):
    result=list(values)
    length=len(result)
    d=1
    while d<length:
        for start in range(0,length,2*d):
            for i in range(start,start+d):
                a,b=result[i],result[i+d]
                result[i],result[i+d]=a+b,a-b
        d*=2
    return result


def basis(vectors):
    pivots={}
    for x in vectors:
        value=x
        while value:
            j=value.bit_length()
            if j not in pivots:
                pivots[j]=value
                break
            value ^= pivots[j]
    return len(pivots)


def check(a,b):
    q=1<<b
    size=1<<(a+b+1)
    count=q**(a+1)
    table=[]
    rows=[]
    bands={}
    for u in range(count):
        coefficients=[u>>(i*b)&(q-1) for i in range(a+1)]
        row=[]
        for x in range(size):
            v=x>>(a+b)
            z=(x>>a)&(q-1)
            y=x% (1<<a)
            center=0
            for i,c in enumerate(coefficients[1:]):
                center |= ((z&c).bit_count()%2)<<i
            sign=v^((z&coefficients[0]).bit_count()%2)^(y==center)
            row.append(1-2*sign)
        table.extend(row)
        spectrum=walsh(row)
        assert sum(map(abs,spectrum))*2**a==size*(3*2**a-4)
        labels=tuple(j for j,w in enumerate(spectrum) if w)
        assert len(labels)==2**a
        rows.append(labels)
        for j in labels:
            bands.setdefault(j,set()).add(u)
    atom_count=Counter(rows)
    assert len(atom_count)==count
    checks=0
    for u,labels in enumerate(rows):
        for candidate in combinations(labels,a+1):
            if basis(x+size for x in candidate)==a+1:
                assert set.intersection(*(bands[x] for x in candidate))=={u}
                checks+=1
    support=[i for i,v in enumerate(walsh(table)) if v]
    s=len(support)
    d=basis(support)
    assert all(len(v)==2**(a*b) for v in bands.values())
    differences=[]
    group={}
    for x in support:
        beta=x%size
        if beta not in group:
            group[beta]=x
        differences.append((x^group[beta])//size)
    t=basis(differences)
    assert (s,t,d)==(2**(a+2*b),b*(a+1),b*(a+1)+a+b+1)
    return dict(a=a,b=b,s=s,t=t,d=d,affine_basis_checks=checks,
                query_count=count,label_repetitions=2**(a*b),
                log_square_equality=s.bit_length()**2==8*t)


if __name__=='__main__':
    print(json.dumps(dict(python=platform.python_version(),
                          cases=[check(a,b) for a,b in [(2,1),(3,1),(2,2),(3,2)]]),indent=2))
