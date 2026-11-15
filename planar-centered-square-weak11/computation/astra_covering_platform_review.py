#!/usr/bin/env python3
"""Independent event-based and exact rectangle checks of covering obstructions."""
from fractions import Fraction as F
from pathlib import Path
from collections import Counter
import argparse, hashlib, json, platform, sys, time

def area_intersection(a,b):
    return max(F(0),min(a[1],b[1])-max(a[0],b[0]))*max(F(0),min(a[3],b[3])-max(a[2],b[2]))

def atom_checks(n):
    L=16*n; h=F(1,L); epsilon=F(1,n)
    source=[s*(F(1,2)+k*h) for s in [-1,1] for k in range(1,L+1)]
    atom_mass=4*h
    assert len(source)==32*n and atom_mass*len(source)==8
    us=sorted({F(j,20) for j in range(-5,6)}|{-h/F(7),h/F(7)})
    vs=[F(-9,5),F(-7,8),F(-5,8),F(0),F(5,8),F(7,8),F(9,5)]
    checked=0
    for u in us:
      for v in vs:
        # True capture radii include both coordinates before sorting; ties counted together.
        events=Counter(max(abs(a-u),abs(v)) for a in source)
        cumulative=F(0)
        for r,multiplicity in sorted(events.items()):
            cumulative+=multiplicity*atom_mass
            assert r>0 and cumulative<=4*r*r,(n,u,v,r,cumulative)
            checked+=1
        assert cumulative==8
    margins=[]
    for u in us:
      if abs(u)<F(1,4):
        old=atom_mass*sum(abs(a-u)<=1 for a in source)
        assert old>=4-4*h
        assert old+epsilon>4
        margins.append(str(old+epsilon-4))
    strong_bound=1-epsilon if n>=4 else None
    return dict(n=n,atoms=len(source),old_mass='8',epsilon=str(epsilon),centers=len(us)*len(vs),capture_events=checked,mixed_margin_min=str(min(F(x) for x in margins)),strong_mixed_area_lower=str(strong_bound) if strong_bound is not None else None)

def fixed_l1_checks():
    supports=[(F(1,2),F(3,2),F(-1,4),F(1,4)),(F(-3,2),F(-1,2),F(-1,4),F(1,4))]
    gsupport=(F(-1,4),F(1,4),F(-1,4),F(1,4))
    epsilons=[F(1,4),F(1,1024)]
    assert sum(8*area_intersection(s,s) for s in supports)==8
    tested=0
    for epsilon in epsilons:
      assert 4*epsilon<=1
      assert 4*epsilon*area_intersection(gsupport,gsupport)==epsilon
      for u in [F(-1,4),F(-1,8),F(0),F(1,8),F(1,4)]:
       for v in [F(-3,4),F(-5,8),F(0),F(5,8),F(3,4)]:
        q=(u-1,u+1,v-1,v+1)
        old=sum(8*area_intersection(q,s) for s in supports)
        new=4*epsilon*area_intersection(q,gsupport)
        assert old==4 and new==epsilon and old+new>4
        tested+=1
    return dict(old_mass='8',epsilons=list(map(str,epsilons)),mixed_witnesses=tested,full_plateau_rectangle_area='3/4',source_free_subregion_area='1/4')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True); args=ap.parse_args()
    start=time.perf_counter()
    cases=[atom_checks(n) for n in [1,3,17]]
    result=dict(status='PASS',python=sys.version,platform=platform.platform(),arithmetic='exact fractions.Fraction; deterministic; no tolerance or random seed',cases=cases,fixed_l1=fixed_l1_checks(),elapsed_seconds=time.perf_counter()-start,scope='Actual finite-event tests cover every radius only at the listed centers. General all-n/all-center proof is analytic. Fixed L1 tests verify exact rectangle integrals, not the global analytic exclusion.')
    root=Path(__file__).resolve().parents[1]
    paths=[root/'proofs/cross-scale-insertion-obstruction.md',root/'proofs/threshold-plateau-insertion.md',root/'proofs/finite-atomic-threshold-platform-insertion.md',Path(__file__).resolve()]
    result['source_sha256']={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
