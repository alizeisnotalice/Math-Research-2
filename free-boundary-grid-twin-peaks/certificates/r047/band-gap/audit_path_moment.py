from collections import deque
from fractions import Fraction


def adjacent(v, n):
    x,y=v
    return tuple((u,w) for u,w in ((x-1,y),(x+1,y),(x,y-1),(x,y+1)) if 0<=u<n and 0<=w<n)


def d1(a,b):
    return abs(a[0]-b[0])+abs(a[1]-b[1])


def growth(n, prefix, excluded=(), limit=None):
    out=list(prefix)
    if len(out)==limit:
        return out
    queue=deque(prefix)
    seen=set(prefix)|set(excluded)
    while queue:
        v=queue.popleft()
        for w in adjacent(v,n):
            if w not in seen:
                out.append(w)
                queue.append(w)
                seen.add(w)
                if limit is not None and len(out)==limit:
                    return out
    return out


def exact_peaks(values,n):
    return {v for v,x in values.items() if all(values[w]<x for w in adjacent(v,n))}


def route(p,q):
    x,y=p
    out=[p]
    while x!=q[0]:
        x += (x<q[0])-(x>q[0])
        out.append((x,y))
    while y!=q[1]:
        y += (y<q[1])-(y>q[1])
        out.append((x,y))
    return out


def audit_case(n,p,q,qrank):
    N=n*n
    path=route(p,q)
    d=len(path)-1
    domain=set(path[:-1])
    prefix=growth(n,[p],set(adjacent(q,n))|{q},qrank-1)
    assert len(prefix)==qrank-1
    original=growth(n,prefix+[q])
    assert len(original)==N and original[qrank-1]==q
    values={v:Fraction(N+1-i,N+1) for i,v in enumerate(original,1)}
    assert exact_peaks(values,n)=={p,q}
    delta=Fraction(1,200)
    begin=max(values[q],1-delta/2)
    end=(begin+values[p])/2
    jcount=2*d+3
    step=(end-begin)/jcount
    assert step>0
    shell={w for v in domain for w in adjacent(v,n)}-domain
    assert len(shell)<=2*d+2
    holes=[(begin+i*step, begin+(i+1)*step) for i in range(jcount)
           if all(not begin+i*step<values[w]<begin+(i+1)*step for w in shell)]
    assert holes
    tests=0
    for low, high in holes:
        seeds={w for w in domain if values[w]>high}
        assert p in seeds and path[-2] not in seeds
        queue=deque(sorted(seeds))
        seen=set(seeds)
        repair=[]
        while queue:
            v=queue.popleft()
            for w in adjacent(v,n):
                if w in domain and w not in seen:
                    seen.add(w)
                    queue.append(w)
                    repair.append(w)
        assert seen==domain
        output=values.copy()
        for i,w in enumerate(repair,1):
            output[w]=high-(high-low)*Fraction(i,len(repair)+1)
        assert exact_peaks(output,n)=={p}
        left = p[0] < n/2
        good=[]
        for v,x in values.items():
            if (v[0]>=n-n//4 if left else v[0]<n//4) and (v[0]+2*v[1])%5==0:
                highest=max(values[w] for w in adjacent(v,n))
                if x<high and highest<=1-delta:
                    assert d1(p,v)>n/8 and d1(q,v)>2
                    good.append((v,highest))
        assert len(good)>=N/40
        v,highest=good[len(good)//2]
        output[v]=(high+highest)/2
        assert exact_peaks(output,n)=={p,v}
        assert max(output,key=output.get)==p
        tests+=1
    return tests


if __name__=="__main__":
    tests=0
    cases=[
        (100,(0,0),(1,1),2),
        (101,(4,6),(7,6),2),
        (100,(54,3),(51,3),120),
        (103,(52,45),(55,47),240),
        (108,(22,50),(28,50),600),
        (110,(63,45),(63,42),2),
        (112,(111,111),(110,109),900),
    ]
    for args in cases:
        tests+=audit_case(*args)
    print(f"exact path repair checks: {len(cases)} input words, {tests} empty windows")
