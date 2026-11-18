"""Finite exact audit of the max-index estimate, singleton constants, and reciprocal symbol."""
from fractions import Fraction as F
from math import comb, factorial
import json,platform

def splits(n,prefix,slots):
    if slots==1:
        yield (*prefix,n)
        return
    for j in range(n+1):
        yield from splits(n-j,(*prefix,j),slots-1)

g0inverse=F(1,2)
V0=F(6)
c=g0inverse*V0
def G(j):
    return g0inverse*(-2)**j/factorial(j)
def h(j):
    return V0*(-1)**j*(j+1)
seen=set()
count=0
for i in range(4):
    p=2*i+1
    for n in range(1,13):
        singles=F(0)
        for js in splits(n,(),p):
            star=max(range(p),key=lambda v:js[v])
            J=js[star]
            volume=F(1)
            for j in js:volume*=factorial(j)
            volume/=factorial(n)
            assert volume<=F(1,2**(n-J))
            if sum(j>0 for j in js)>=2:
                seen.add(('max_G' if star<i else 'max_V','tie' if js.count(J)>1 else 'unique'))
            if sum(j>0 for j in js)==1:
                term=volume
                for l,j in enumerate(js):term*=G(j) if l<i else h(j)
                singles+=term
            count+=1
        assert singles==(i+1)*c**i*h(n)+(i*g0inverse**(i-1)*V0**(i+1)*G(n) if i else 0)
def S(k):
    return sum((c/2)**i*h(k-2*i)/factorial(i) for i in range(k//2+1))
for k in range(31):
    assert h(k)==sum((-c/2)**i*S(k-2*i)/factorial(i) for i in range(k//2+1))
print(json.dumps({'python':platform.python_version(),'exact_max_index_compositions':count,
  'mixed_classes_seen':sorted([list(p) for p in seen]),'singleton_checks':48,
  'exact_symbol_inverse_orders':31},indent=2))
