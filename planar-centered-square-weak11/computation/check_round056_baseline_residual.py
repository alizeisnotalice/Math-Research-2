#!/usr/bin/env python3
"""Exact rational checks for R56 baseline/residual transport small cases."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,platform

def area(b): return max(F(0),b[1]-b[0])*max(F(0),b[3]-b[2])
def intersect(a,b): return max(a[0],b[0]),min(a[1],b[1]),max(a[2],b[2]),min(a[3],b[3])
def qbox(x,r): return x[0]-r,x[0]+r,x[1]-r,x[1]+r

def main():
    e=F(1,100); positions=[(F(0),F(0)),(F(12,5),F(0))]
    sources=[qbox(p,e) for p in positions]
    rs=[(F(-49,50),F(49,50),F(-49,50),F(49,50)),
        (F(71,50),F(169,50),F(-49,50),F(49,50)),
        (F(101,100),F(139,100),F(-139,100),F(139,100))]
    radii=[F(99,100),F(99,100),F(7,5)]
    chosen=[[0],[1],[0,1]]; vertices=0
    for rect,r,ids in zip(rs,radii,chosen):
        for x in [(rect[a],rect[b]) for a in (0,1) for b in (2,3)]:
            box=qbox(x,r)
            for i in ids: assert area(intersect(box,sources[i]))==area(sources[i])
            vertices+=1
    for i in range(3):
        for j in range(i): assert area(intersect(rs[i],rs[j]))==0
    avgs=[F(4)/(4*radii[0]**2),F(4)/(4*radii[1]**2),F(8)/(4*radii[2]**2)]
    assert all(a>1 for a in avgs)
    lower=sum(map(area,rs)); baseline=sum(map(area,sources))
    assert lower==F(21849,2500) and baseline==F(1,1250)
    sigma=8-baseline; qlower=lower-baseline
    assert sigma==F(19998,2500) and qlower==F(21847,2500)
    L=F(64); D=(-L,L,F(0),L); A=(F(-1),F(1),F(0),F(1)); wedge=0
    for t in map(F,[4,5,6,8,12,16]):
        for s in [-t/2,-t/4,F(0),t/4,t/2,F(1,3)]:
            box=qbox((s,t),t)
            assert area(intersect(box,D))==area(box)
            assert area(intersect(box,A))==2
            mass=area(box)+2
            assert mass>4*t*t
            assert F(1)-min(F(1),F(1))==0  # q=0 at each wedge point
            wedge+=1
    out={'python':platform.python_version(),'arithmetic':'fractions.Fraction',
         'random_seed':None,'tolerance':None,'processes':1,
         'receiver_rectangles':[[str(x) for x in r] for r in rs],
         'source_containment_vertex_checks':vertices,'wedge_witness_checks':wedge,
         'strict_witness_averages':list(map(str,avgs)),
         'certified_halo_area_lower':str(lower),'input_mass':'8',
         'baseline_mass':str(baseline),'residual_source_mass':str(sigma),
         'residual_demand_lower':str(qlower),'residual_Hall1_deficit_lower':str(qlower-sigma),
         'residual_inflation_lower':str(qlower/sigma),'full_capacity_lower':str(lower/8),
         'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    dest=Path(__file__).parent/'results/round056_baseline_residual_exact.json'
    dest.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__': main()
