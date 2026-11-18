"""All linear quotient kernels in dimensions 3,5,7; integer signed spectra."""
from collections import Counter
from itertools import combinations
from multiprocessing import Pool
import json
import platform


def images(basis):
    out=[(0,0)]
    for i, v in enumerate(basis):
        out += [(a^v, tag|(1<<i)) for a,tag in out]
    return dict(out)


def kernels(n):
    for k in range(n+1):
        for piv in combinations(range(n), k):
            free=[(i,j) for i,p in enumerate(piv) for j in range(p) if j not in piv]
            for mask in range(1<<len(free)):
                rows=[1<<p for p in piv]
                for k0,(i,j) in enumerate(free):
                    if (mask>>k0)&1: rows[i]|=1<<j
                yield rows


def rank(vs):
    pivots={}
    for v in vs:
        while v and v.bit_length() in pivots:
            v ^= pivots[v.bit_length()]
        if v:
            pivots[v.bit_length()] = v
    return len(pivots)


def walsh(row):
    row=row.copy()
    w=1
    while w<len(row):
        for i in range(0,len(row),2*w):
            for j in range(i,i+w):
                a,b=row[j],row[j+w]
                row[j],row[j+w]=a+b,a-b
        w*=2
    return row


def checks(task):
    rows, coeff, P, Q, e = task
    coords=images(rows)
    def split(a):
        rep=a
        for v in reversed(rows):
            if (rep>>(v.bit_length()-1))&1: rep ^= v
        return rep, coords[a^rep]
    u,v=({split(a)[0] for a in X} for X in (P,Q))
    A,B=sorted((len(u),len(v)),reverse=True)
    bands={}
    for a,c in coeff.items():
        rep,k=split(a)
        bands.setdefault(rep,{})[k]=c
    M=0
    depth=[0]*len(coords)
    pure=[[],[]]
    regular=set()
    for rep, dic in bands.items():
        base=next(iter(dic))
        fill=1 << rank(x^base for x in dic)
        M+=fill
        activity={i for i,c in enumerate(walsh([dic.get(i,0) for i in range(len(coords))])) if c}
        for i in activity:depth[i]+=1
        if rep not in (0,split(e)[0]):
            regular |= activity
        if rep not in u or rep not in v:
            pure[int(rep not in u)].append(activity)
            assert len(activity)*fill==len(coords)
    assert all(all(x==ss[0] for x in ss) for ss in pure)
    N=(len(P)//2)
    assert len(regular)*2*N <= (A+B)*len(coords)
    return M*max(depth)


def main():
    summary=[]
    for m in (1,2,3):
        N=1<<m
        n=2*m+1
        e=1<<(2*m)
        P={x|(t*e) for t in range(2) for x in range(N)}
        Q={(y<<m)|(t*e) for t in range(2) for y in range(N)}
        row=[-1 if (z & (N-1) if not z&e else (z>>m)&(N-1))==0 else 1 for z in range(1<<n)]
        coeff={i:c for i,c in enumerate(walsh(row)) if c}
        expected=(P|Q)-{e}
        if m==1:expected.remove(0)
        assert set(coeff)==expected
        with Pool(12) as pool:
            products=pool.map(checks, ((h,coeff,P,Q,e) for h in kernels(n)),chunksize=100)
        summary.append(dict(m=m,N=N,s=len(coeff),d=rank(coeff),kernels=len(products),
                            minimum_actual_LM=min(products),max_actual_LM=max(products)))
    print(json.dumps(dict(python=platform.python_version(),cases=summary),indent=2))


if __name__=='__main__':
    main()
