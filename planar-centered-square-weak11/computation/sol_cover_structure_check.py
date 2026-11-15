from fractions import Fraction as F
import hashlib,json,platform
from pathlib import Path

def area(rects):
    xs=sorted({x for r in rects for x in r[:2]}); out=F(0)
    for l,r in zip(xs,xs[1:]):
        ys=sorted((b,t) for a,z,b,t in rects if a<r and z>l)
        total=F(0); end=None
        for b,t in ys:
            if end is None or b>end: total+=t-b;end=t
            elif t>end: total+=t-end;end=t
        out+=(r-l)*total
    return out
cases=[]
for n in (2,8,32):
    d=F(1,4*n); e=8*d+4*d*d
    rects=[]
    for i in range(n):
        for j in range(n):
            x=F(2*i+1,2*n)-F(1,2);y=F(2*j+1,2*n)-F(1,2)
            rects.append((x-d/2,x+d/2,y-d/2,y+d/2))
    new=[(a-d,b+d,c-d,z+d) for a,b,c,z in rects]
    oldarea=area(rects);newarea=area(new)
    assert oldarea==F(1,16) and newarea==F(9,16)
    assert (1+d)**2==1+e/4
    cases.append(dict(n=n,epsilon=str(e),delta=str(d),old=str(oldarea),new=str(newarea),gain=str(newarea-oldarea),ratio=str((newarea-oldarea)/e)))
root=Path(__file__).resolve().parents[1]
names=['proofs/general-dual-interface.md','proofs/kernel-family-floor.md','proofs/cross-scale-insertion-obstruction.md','proofs/threshold-plateau-insertion.md']
print(json.dumps(dict(python=platform.python_version(),algorithm='exact Fraction sweep union; deterministic single process; no random seed or tolerance',cases=cases,source_sha256={s:hashlib.sha256((root/s).read_bytes()).hexdigest() for s in names}),indent=2))
