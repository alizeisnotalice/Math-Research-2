#!/usr/bin/env python3
"""Small deterministic floating rectangle-union probes; no global claims."""
from pathlib import Path
from itertools import combinations, product
from fractions import Fraction as F
import argparse, hashlib, json, math, platform, random, sys, time
TOPIC=Path(__file__).resolve().parents[1]

def rectangles(atoms):
    out=[]
    for mask in range(1,1<<len(atoms)):
        group=[a for j,a in enumerate(atoms) if mask>>j&1]
        mass=sum(a[2] for a in group)
        s=math.sqrt(float(mass))/2
        x0=max(float(a[0]) for a in group)-s
        x1=min(float(a[0]) for a in group)+s
        y0=max(float(a[1]) for a in group)-s
        y1=min(float(a[1]) for a in group)+s
        if x0<x1 and y0<y1:
            out.append((x0,x1,y0,y1))
    return tuple(set(out))

def area(rects):
    xs=sorted(set(t for a,b,c,d in rects for t in (a,b)))
    total=0.0
    for a,b in zip(xs,xs[1:]):
        mid=(a+b)/2
        ints=sorted((c,d) for x0,x1,c,d in rects if x0<mid<x1)
        length=0.0; lo=None; hi=None
        for c,d in ints:
            if lo is None: lo,hi=c,d
            elif c>hi: length+=hi-lo;lo,hi=c,d
            else: hi=max(hi,d)
        if lo is not None: length+=hi-lo
        total+=(b-a)*length
    return total

def trial(old,eps,kind):
    nu=((F(0),F(0),eps),)
    r_old=rectangles(old);r_nu=rectangles(nu);r_sum=rectangles(old+nu)
    m=sum(a[2] for a in old)
    hsum=area(r_sum);hunion=area(r_old+r_nu)
    pure=max(0.0,hsum-hunion)
    return {'kind':kind,'old':[[str(t) for t in a] for a in old],
            'epsilon':str(eps),'M':str(m),'pure_mixed_area_float':pure,
            'pure_over_M_float':pure/float(m),
            'gap_over_M_float':pure-float(m),
            'gap_over_M_plus_4eps_float':pure-float(m+4*eps)}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();start=time.perf_counter()
    # Exact root comparisons for the two-old-atom halo counterexample.
    assert F(2)>F(7,5)**2 and F(7,5)>0
    two_old=((F(-1,10),F(0),F(1)),(F(1,10),F(0),F(1)))
    exact_halo_float=2-math.sqrt(2)/5
    assert abs(area(rectangles(two_old))-exact_halo_float)<1e-12
    menu=[]
    seed=20261115
    rng=random.Random(seed)
    epses=(F(1,100),F(1,10),F(1),F(10))
    for index in range(100):
        count=2+index%4
        coords=set()
        while len(coords)<count:
            q=(F(rng.randrange(-8,9),4),F(rng.randrange(-8,9),4))
            if q!=(F(0),F(0)):
                coords.add(q)
        coords=tuple(sorted(coords))
        raw=[rng.randrange(1,9) for _ in coords]
        denom=sum(raw)
        old=tuple((q[0],q[1],F(w,denom)) for q,w in zip(coords,raw))
        eps=epses[(index//4)%4]
        row=trial(old,eps,'seeded-small-input')
        row['index']=index
        menu.append(row)
    best=max(menu,key=lambda c:c['pure_over_M_float'])
    violations=[c for c in menu if c['gap_over_M_float']>1e-10]
    output={'status':'completed','python':sys.version,'platform':platform.platform(),
            'scope':'deterministic menu only; floating union-area exploration, not a global or rigorous negative certificate',
            'arithmetic':'rational input, IEEE double square roots and union areas; exploratory detection tolerance 1e-10',
            'source_counts':'2, 3, 4, or 5 old atoms plus the new point atom at origin',
            'case_count':len(menu),'seed':seed,'epsilon_over_M_menu':['1/100','1/10','1','10'],
            'coordinate_menu':'{-8,...,8}/4 in each coordinate, distinct nonzero old positions; positive integer weights 1..8 normalized to M=1',
            'best_pure_over_M':best,'detected_M_violations':violations,'cases':menu,
            'two_atom_halo_counterexample':{'M':'2','area_exact':'2-sqrt(2)/5','area_float':exact_halo_float,
              'root_comparison_exact':'2>49/25, hence sqrt(2)>7/5'},
            'source_sha256':{str(p.relative_to(TOPIC)):hashlib.sha256(p.read_bytes()).hexdigest() for p in
               (TOPIC/'proofs/joint-old-reserve-counterexample-notes.md',Path(__file__).resolve())},
            'elapsed_seconds':time.perf_counter()-start}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output,indent=2))
if __name__=='__main__': main()
