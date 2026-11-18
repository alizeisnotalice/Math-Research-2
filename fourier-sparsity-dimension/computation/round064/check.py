"""Exact test of all 4-variable truth tables with a fixed 2-dimensional core."""
from collections import Counter
from fractions import Fraction
import itertools
import json
import platform

def fwht(a):
    a=list(a)
    w=1
    while w<len(a):
        for i in range(0,len(a),2*w):
            for j in range(i,i+w):
                a[j],a[j+w]=a[j]+a[j+w],a[j]-a[j+w]
        w*=2
    return a

def rank(xs):
    p={}
    for v in xs:
        while v.bit_length() in p:
            v^=p[v.bit_length()]
        if v:
            p[v.bit_length()]=v
    return len(p)

def check():
    physical=4
    windows={}
    for k in range(3):
        samples=Counter()
        for b,*vec in itertools.product(range(physical),repeat=k+1):
            span={b}
            for v in vec:
                span|={x^v for x in span}
            samples[sum(1<<x for x in span)]+=1
        windows[k]=samples
    row=[fwht([1-2*((mask>>y)&1) for y in range(4)]) for mask in range(16)]
    rank_of=[rank(i for i in range(4) if (s>>i)&1) for s in range(16)]
    stats=Counter()
    comparisons=0
    for f in range(65536):
        rows=[row[(f>>(4*x))&15] for x in range(4)]
        coeff=[fwht(v) for v in zip(*rows)]
        ns=[sum(a!=0 for a in v) for v in coeff]
        active=[sum(1<<x for x in range(4) if rows[x][b]) for b in range(4)]
        pattern=tuple(active),tuple(ns)
        stats[pattern]+=1
    for (active,ns),multiplicity in stats.items():
        for k,w in windows.items():
            M=1<<k
            total=sum(w.values())
            survival=[]
            available=sum(1<<b for b,u in enumerate(ns) if u)
            r=rank_of[available]
            numer=0
            for sample,count in w.items():
                alive=sum(1<<b for b in range(4) if active[b]&sample)
                numer+=rank_of[alive]*count
            for b,n in enumerate(ns):
                if n:
                    p=Fraction(active[b].bit_count(),4)
                    loss=Fraction(sum(v for key,v in w.items() if not(active[b]&key)),total)
                    e=min(Fraction(n-1,2*M),(1-p)/(1+(M-1)*p),(1-p)**(k+1))
                    assert loss<=e
                    survival.append(1-e)
                else:
                    survival.append(0)
            best=max((sum(survival[b] for b in range(4) if I>>b&1)
                      for I in range(16) if I&available==I and I.bit_count()==rank_of[I]==r),default=0)
            assert numer>=best*total
            comparisons+=multiplicity
    return {"python":platform.python_version(),"truth_tables":65536,"ambient":4,"kernel_dim":2,
            "M":[1,2,4],"affine_maps_per_M":[sum(s.values()) for s in windows.values()],
            "distinct_fiber_patterns":len(stats),"table_by_scale_cases":comparisons,
            "exact_real_rational":True,"scope":"all band-loss and expected projected basis-rank inequalities"}
if __name__=="__main__":
    print(json.dumps(check(),indent=2))
