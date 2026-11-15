"""Independent Burnside/quotient-graph audit of D4 independent-set layers."""
from functools import cache
from pathlib import Path
import json


def raw_layers(n):
    valid=[s for s in range(1<<n) if not s&(s<<1)];f={(0,0):1}
    for _ in range(n):
        g={}
        for (previous,k),v in f.items():
            for s in valid:
                if not s&previous:
                    key=(s,k+s.bit_count());g[key]=g.get(key,0)+v
        f=g
    out=[0]*(n*n+1)
    for (_,k),v in f.items():out[k]+=v
    return tuple(out)


def fixed_layers(n,op):
    N=n*n
    def transform(v):
        x,y=v%n,v//n
        if op>=4:x=n-1-x
        for _ in range(op%4):x,y=n-1-y,x
        return y*n+x
    cycles=[];which={}
    for v in range(N):
        if v in which:continue
        orbit=[];w=v
        while w not in which:
            which[w]=len(cycles);orbit.append(w);w=transform(w)
        cycles.append(orbit)
    qadj=[0]*len(cycles);invalid=0
    for v in range(N):
        for w in (v+1,v+n):
            if w>=N or abs(v%n-w%n)+abs(v//n-w//n)!=1:continue
            a,b=which[v],which[w]
            if a==b:invalid|=1<<a
            else:qadj[a]|=1<<b;qadj[b]|=1<<a
    @cache
    def poly(mask):
        if not mask:return (1,)+(0,)*N
        bit=mask&-mask;v=bit.bit_length()-1;weight=len(cycles[v])
        exclude=poly(mask^bit);include=poly((mask^bit)&~qadj[v]);out=list(exclude)
        for k in range(N+1-weight):out[k+weight]+=include[k]
        return tuple(out)
    return poly(((1<<len(cycles))-1)^invalid)

n=7;raw=raw_layers(n);fixed=[raw]+[fixed_layers(n,op) for op in range(1,8)]
orbits=[]
for k in range((n*n+1)//2+1):
    total=sum(p[k] for p in fixed);q,rem=divmod(total,8);assert rem==0;orbits.append(q)
z=json.loads((Path(__file__).parent/'results/n7-exact-r005.json').read_text())
assert orbits==z['representatives_by_size']
assert sum(raw)==z['raw_independent_sets']==1280128950
print(json.dumps({'verified':True,'n':n,'raw_layers':raw,'fixed_sets_by_symmetry':[sum(p) for p in fixed],
                  'D4_orbits_by_size':orbits,'D4_orbits_total':sum(orbits)},indent=2))
