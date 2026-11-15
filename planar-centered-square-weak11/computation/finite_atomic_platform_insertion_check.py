#!/usr/bin/env python3
"""Exact, single-process checks of the finite atomic critical platform."""
from fractions import Fraction as F
from pathlib import Path
import argparse, hashlib, json, platform, sys, time

TOPIC = Path(__file__).resolve().parents[1]

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def check(n):
    eps=F(1,n); L=16*n; h=F(1,L)
    positions=tuple(s*(F(1,2)+k*h) for s in (-1,1) for k in range(1,L+1))
    mass=4*h
    assert len(set(positions))==32*n and mass>0 and mass*len(positions)==8
    # Rational centers include strip edges and capture-event phases.
    us=tuple(sorted(set((F(-1,4),F(-1,8),-h/3,F(0),h/3,F(1,8),F(1,4)))))
    vs=(F(0),F(1,2),F(5,8),F(-5,8),F(3,4),F(3,2))
    event_total=0; old_center_total=0; mixed_center_total=0
    for u in us:
        assert abs(u)<=F(1,4)
        for v in vs:
            # A step in captured mass can occur only at one of these radii.
            events=sorted(set([abs(v)]+[abs(t-u) for t in positions]))
            for r in events:
                if r<=0:
                    continue
                old=mass*sum(abs(t-u)<=r for t in positions) if r>=abs(v) else F(0)
                assert old<=4*r*r, (n,u,v,r,old)
                event_total+=1
            old_center_total+=1
        # Fixed r=1 and v=±5/8 are interior mixed witnesses for |u|<1/4.
        if abs(u)<F(1,4):
            cplus=int((F(1,2)+u)//h)
            cminus=int((F(1,2)-u)//h)
            assert 0<cplus<L and 0<cminus<L
            assert cplus+cminus in (L-1,L)
            old=mass*(cplus+cminus)
            assert old>=4-8*h
            assert old+eps>=4+eps/2>4
            for v in (F(-5,8),F(5,8)):
                assert abs(u)<1 and abs(v)<1
                assert 4*v*v>eps  # no single-atom strict halo at this center
                mixed_center_total+=1
    # Three exact polynomial identities and endpoint checks, not symbolic proofs.
    for r in (F(0),F(1,4),F(1,2),F(3,4),F(1),F(3,2),F(4)):
        assert r*r-r+F(1,4)==(r-F(1,2))**2
        assert r*r-2*r+1==(r-1)**2
        assert 4*r*r-4*(r-F(1,4))==4*(r-F(1,2))**2
    area=F(1,2)*2*F(1,4)
    assert area==F(1,4)
    return {'n':n,'epsilon':str(eps),'L':L,'h':str(h),'old_atoms':len(positions),
            'old_atom_mass':str(mass),'old_total_mass':'8',
            'old_all_radius_centers':old_center_total,'old_capture_events_checked':event_total,
            'mixed_centers_checked':mixed_center_total,'mixed_margin_at_least':str(eps/2),
            'pure_mixed_area_lower':'1/4','status':'PASS'}

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    start=time.perf_counter()
    cases=[check(n) for n in (1,2,7,32)]
    source_paths=[TOPIC/'proofs/finite-atomic-threshold-platform-insertion.md',
                  TOPIC/'proofs/cross-scale-insertion-obstruction.md',
                  TOPIC/'proofs/threshold-plateau-insertion.md', Path(__file__).resolve()]
    result={'status':'PASS','python':sys.version,'platform':platform.platform(),
            'arithmetic':'fractions.Fraction exact rational; no random inputs or floating tolerances',
            'cases':cases,'total_event_checks':sum(c['old_capture_events_checked'] for c in cases),
            'source_sha256':{str(p.relative_to(TOPIC)):sha(p) for p in source_paths},
            'elapsed_seconds':time.perf_counter()-start,
            'scope':'Finite deterministic source-distance events validate all radii at each listed center; all-center/all-n theorem is analytic.'}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False,indent=2))

if __name__=='__main__':
    main()
