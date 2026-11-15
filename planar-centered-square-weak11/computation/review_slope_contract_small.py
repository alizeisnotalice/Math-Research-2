"""Bounded exact checks for the local-slope contract; no true-path geometry run."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from datetime import datetime,timezone
import argparse,hashlib,json,sys,time

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();start=time.perf_counter();cases=points=0;minmargin=None
    for lam,d0,H,ratio,u in product((F(1),F(3,2)),(F(1),F(3),F(6)),(F(1,100),F(1,4)),(F(1,2),F(99,100)),(-1,1)):
        m=lam*d0*d0*ratio;q=lam*d0*d0-m;L=lam*abs(u)*(2*abs(d0)+abs(u)*H);delta=min(H,q/(2*L));assert delta>0;cases+=1
        for j in range(-16,17):
            h=delta*F(j,16);d=d0+u*h;residual=lam*d*d-m
            assert d>0 and abs(lam*d*d-lam*d0*d0)<=L*abs(h)and residual>=q/2
            minmargin=residual-q/2 if minmargin is None else min(minmargin,residual-q/2);points+=1
    # A legitimate R=3 cap contact: full-box radius 3, x span 6+s.
    # At s=0 it is empty; for negative s it has positive x width.
    cap_width=lambda s:max(F(0),6-(6+s))
    assert cap_width(F(0))==0 and cap_width(F(-1,100))==F(1,100)and cap_width(F(1,100))==0
    # Fixed-y clipping can make a bbox permanently empty even when its span is 0.
    # Redundant periodic source at y=13, radius1, receiving base y-cell [0,10].
    assert max(F(0),min(F(10),F(14))-max(F(0),F(12)))==0
    # A raw left endpoint can cross F_left even though its earlier clipped value
    # is the persistent F boundary. Its velocity changes from 0 to 1 at s=0.
    clipped_left=lambda s:max(F(0),s)
    assert clipped_left(F(-1,4))==clipped_left(F(-1,8))==0 and clipped_left(F(1,8))==F(1,8)
    out={'status':'passed','scope':'finite exact margin and edge-case checks only, not a real 204-label slope certificate','rational_empty_margin_cases':cases,'rational_point_checks':points,'minimum_residual_above_q_over_2':str(minmargin),'span_2R_center_wall_detected':True,'fixed_y_clipped_empty_detected':True,'unselected_raw_clip_branch_event_detected':True,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'python':sys.version,'command':[sys.executable,*sys.argv],'verified_utc':datetime.now(timezone.utc).isoformat(),'elapsed_seconds':time.perf_counter()-start}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
