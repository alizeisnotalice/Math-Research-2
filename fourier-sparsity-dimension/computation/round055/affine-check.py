"""All n=4 truth tables, rank r evaluators in two affine-coordinate normal forms."""
from collections import Counter
from fractions import Fraction
import json
import platform

def rank(xs):
    p = {}
    for x in xs:
        while x.bit_length() in p:
            x ^= p[x.bit_length()]
        if x:
            p[x.bit_length()] = x
    return len(p)

n = 4
N = 1 << n
masks = [sum(1 << x for x in range(N) if (x&a).bit_count() % 2) for a in range(N)]
forms = []
for r in range(1,n+1):
    for constant in (False, True):
        d = r-int(constant)
        labels = [(a & ((1<<d)-1)) + (int(constant)<<d) for a in range(N)]
        forms.append((r,constant,labels))
counts = {(r,c):Counter() for r,c,_ in forms}
for f in range(1<<N):
    w = [N-2*(f^m).bit_count() for m in masks]
    S = [a for a in range(N) if w[a]]
    for r,constant,labels in forms:
        if rank(labels[a] for a in S)!=r:
            continue
        low = [a for a in S if labels[a]]
        deficit = sum(abs(w[a]) for a in low)-2*N+(2*N>>r)
        assert deficit >= 0
        counts[r,constant][deficit]+=1
        if deficit:
            least = {1:8,2:8,3:4,4:10}[r]
            assert deficit>=least
            continue
        if constant:
            assert len(S)==r==1
            continue
        v0 = low[0]
        a = abs(w[v0])
        assert all(abs(w[x])==a for x in low)
        assert any(all(((w[x]>0)==(w[v0]>0))==(((x^v0)&t).bit_count()%2==0)
                       for x in low) for t in range(N))
        cells = {k:[x for x in low if labels[x]==k] for k in range(1,1<<r)}
        anchor=cells[1][0]
        K={x^anchor for x in cells[1]}
        assert len(K)==1<<rank(K)
        assert all({x^cell[0] for x in cell}==K for cell in cells.values())
        if r>=2:
            base=cells[1][0]^cells[2][0]^cells[3][0]
            deleted={x^base for x in K}
            assert len(K)*(1<<r)==1<<rank([x^base for x in low])
            assert max(abs(w[x]) for x in deleted)>=((1<<(r-1))-1)*a
out=[]
for (r,c),hist in counts.items():
    out.append(dict(r=r,contains_constant=c,faithful=sum(hist.values()),endpoint=hist[0],
                  positive_gap=str(Fraction(min(x for x in hist if x),N)) if any(hist) else None))
print(json.dumps(dict(python=platform.python_version(),ambient_dimension=n,truth_tables=1<<N,
                      models=out),indent=2))
