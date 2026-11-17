"""Small exact root-screen audit from the independently rebuilt 67 S1 mothers.

No production root-screen code is imported.  Coefficients are reconstructed
via H=(F+a x)(F+(T-a)x)=F^2+T x F+N x^2.  Both Schur scalings are explicit.
"""
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

OWN=Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/algebraic_global/round055-s1-independent/full-v1')
PROD=Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/section_eigen39/round055-degree8-entry/root-and-trees-v1')
POINTS=((1,1),(1,2),(2,1),(1,4),(4,1),(3,4),(4,3),(6,5),(5,6))

def conv(a,b):
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): out[i+j]+=x*y
    return out

def H(row):
    T,N,b,c,d,e,f,g=row
    F=[1,0,b,c,d,e,f,g,1]
    h=conv(F,F)
    for i,x in enumerate(F): h[i+1]+=T*x
    h[2]+=N                         # essential degree-2 norm contribution
    assert len(h)==17 and h[0]==h[-1]==1
    return h

def schur_desc(a):
    """Exact strict Schur descent; input coefficients highest-to-constant."""
    a=list(a)
    while len(a)>1:
        d=math.gcd(*a)
        if a[0]<0: d=-d
        a=[x//d for x in a]
        if abs(a[-1])>=a[0]: return len(a)-1
        lead,last=a[0],a[-1]; deg=len(a)-1
        a=[lead*a[i]-last*a[deg-i] for i in range(deg)]
    return None

def screen(row):
    h=H(row)
    for j,(p,q) in enumerate(POINTS):
        if sum(h[i]*p**i*q**(16-i) for i in range(17))<=0: return 10+j
    # H(2x), largest coefficient first.
    fail=schur_desc([h[i]*(2**i) for i in range(16,-1,-1)])
    if fail is not None: return 30+fail
    # rev(H)(2x), largest coefficient first; this is 2^(16-i), not 2^(14-i).
    fail=schur_desc([h[i]*(2**(16-i)) for i in range(17)])
    if fail is not None: return 50+fail
    return 255

def main():
    rows=set()
    for p in sorted((OWN/'cores').glob('*.json')):
        x=json.loads(p.read_text()); q=x['summary']
        for b,c,d,e,f,g,*_ in x['row7_records']: rows.add((q['T'],q['N'],b,c,d,e,f,g))
    assert len(rows)==67
    expected={tuple(x['binding']['row']):x['screen']['code'] for p in PROD.glob('mother-*.json') for x in [json.loads(p.read_text())]}
    assert len(expected)==67 and set(expected)==rows
    got={row:screen(row) for row in sorted(rows)}
    assert got==expected
    out={'status':'passed_67_exact_root_codes','rows':67,'codes':dict(sorted(Counter(got.values()).items())),
         'own_prefix_summary_sha256':hashlib.sha256((OWN/'summary.json').read_bytes()).hexdigest(),
         'production_mother_payload_sha256':hashlib.sha256(b''.join(p.read_bytes() for p in sorted(PROD.glob('mother-*.json')))).hexdigest(),
         'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (Path(__file__).with_name('root-screen-independent-audit.json')).write_text(json.dumps(out,sort_keys=True,indent=2)+'\n')
    print(json.dumps(out,sort_keys=True))

if __name__=='__main__': main()
