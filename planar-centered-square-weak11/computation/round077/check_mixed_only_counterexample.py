#!/usr/bin/env python3
"""R77 exact radical bounds for the analytic mixed-only counterexample.

No random sample, receiver integration or producer's full-halo replay.
"""
from fractions import Fraction as F
import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys
import time


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();start=time.perf_counter();raw=a.input.read_bytes();d=json.loads(raw)
    m=F(d['internal_mass']);left=F(d['left_distance']);bottom=F(d['bottom_distance']);b1,b2=map(F,d['outside_masses'])
    assert (m,b1,b2,left,bottom)==(1,4,8,F(28,25),F(161,50))
    r1sq=(m+b1)/4;single1sq=b1/4;single2sq=b2/4;Tsq=(m+b1+b2)/4;Rsq=(b1+b2)/4
    Tupper=F(d['triple_radius_upper']);Rlower=F(d['residual_pair_radius_lower']);s2upper=F(d['bottom_singleton_radius_upper'])
    assert 1<r1sq<left*left<single2sq<Rsq<Tsq
    assert Tupper*Tupper>Tsq and Rlower*Rlower<Rsq and s2upper*s2upper>single2sq
    assert Tupper+s2upper<bottom<2*Rlower
    assert bottom>3 and left<2
    A1lower=F(3,2)+left-(1+left)*left
    triplelower=F(1,2)-bottom*(Tupper-Rlower)
    total_lower=A1lower+triplelower
    assert A1lower==F(307,1250)and triplelower==F(13569,50000)
    assert total_lower==F(25849,50000)>m/2
    out={'schema':'r077-mixed-only-counterexample-checks-v1','all_passed':True,
         'sources':{'positions':[['0','0'],[str(-left),'0'],['0',str(-bottom)]],'masses':['1','4','8']},
         'mixed_only_exact_area':'78/25-53sqrt(5)/50-161sqrt(13)/100+161sqrt(3)/50',
         'first_pair_strict_lower':str(A1lower),'triple_strict_lower':str(triplelower),
         'mixed_only_strict_lower':str(total_lower),'strict_excess_above_half':str(total_lower-m/2),
         'geometric_separation':'A2 empty; D>U+sqrt(2); first mixed pair and triple separated; triple meets neither residual singleton',
         'scope':'exact rational inequalities for analytic geometry; no complete rectangle union, random search or global-grid qualification performed here',
         'receipt':{'python':sys.version,'platform':platform.platform(),'input_sha256':hashlib.sha256(raw).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'wall_seconds':time.perf_counter()-start}}
    a.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k]for k in ['all_passed','mixed_only_strict_lower','strict_excess_above_half']}))


if __name__=='__main__':main()
