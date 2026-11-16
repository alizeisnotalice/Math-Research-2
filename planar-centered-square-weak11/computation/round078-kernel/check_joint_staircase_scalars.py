#!/usr/bin/env python3
"""R78 actual-prefix source binding and staircase scalar checks.

No halo union, root approximation, position search or producer menu replay.
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
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input',type=Path,required=True);p.add_argument('--actual-input',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();start=time.perf_counter();raw=a.input.read_bytes();d=json.loads(raw)
    actualraw=a.actual_input.read_bytes();actual=json.loads(actualraw)
    assert hashlib.sha256(actualraw).hexdigest()==d['actual_input_sha256']
    n=d['n'];assert n==20;K=n//2;e=F(1,n*n);h=F(n+1,n*n)
    expected=[(i*h,j*h)for i in range(n)for j in range(n)]+[tuple(map(F,p))for p in d['outside_positions']]
    masses=[e]*(n*n)+list(map(F,d['outside_masses']))
    assert [tuple(map(F,p))for p in actual['positions']]==expected
    assert list(map(F,actual['masses']))==masses and actual['S_indices']==list(range(n*n))
    left=-expected[-2][0];bottom=-expected[-1][1]
    assert left==F(71,200)and bottom==F(171,200)
    V2=F(3,16);U1sq=F(301,1600);U9sq=F(309,1600);U10sq=F(310,1600)
    Lbound=F(d['residual_pair_radius_lower']);Hbound=F(d['joint_radius_upper'])
    assert Lbound*Lbound<V2<U1sq<=U9sq<Hbound*Hbound<F(1,4)
    assert left<Lbound and bottom<2*Lbound and bottom-F(1,2)>F(1,4)
    assert F(1,8)<left*left  # residual beta2 radius sqrt(2)/4 < left
    assert Lbound>8*h+e and U10sq<(9*h+e)**2
    assert left+9*h<bottom and F(301,400)>bottom*bottom
    vertical_null_base=(bottom+h)**2-(F(3,4)+F(2,40))
    vertical_null_increment=2*h*(bottom+h)+h*h-F(1,40)
    assert vertical_null_base>0 and vertical_null_increment>0
    assert F(1,1600)/(2*Lbound)<h
    delta_first_lower=F(1,1600)/(2*Hbound)
    first_upper=F(1,800)-(bottom+2*e)*delta_first_lower
    widths_sum_upper=8*(Hbound-e)-36*h
    increment_upper=widths_sum_upper/(1600*Lbound)
    total=first_upper+increment_upper
    assert delta_first_lower==F(1,1408)and first_upper==F(9,14080)
    assert widths_sum_upper==F(161,100)and increment_upper==F(161,68800)
    assert total<F(d['claimed_joint_upper'])
    out={'schema':'r078-joint-staircase-scalar-result-v1','all_passed':True,'actual_source_count':len(expected),'actual_input_matches_formula':True,
         'menu_upper_count':4*K*K+3,'nonempty_old_prefix_subset_labels':4*K*K,'pure_outside_labels':3,
         'proved_active_joint_prefixes':[[k,1]for k in range(1,10)],
         'vertical_null_base_gap':str(vertical_null_base),'vertical_null_increment_lower':str(vertical_null_increment),
         'first_survivor_strict_upper':str(first_upper),'remaining_eight_increments_strict_upper':str(increment_upper),
         'joint_survivor_strict_upper':str(total),'simple_strict_upper':d['claimed_joint_upper'],
         'gap_to_simple_upper':str(F(d['claimed_joint_upper'])-total),
         'scope':'scalar proof checks and exact raw-input binding only; no 403-menu rebuild, root enclosure, halo integration or independent area replay',
         'receipt':{'python':sys.version,'platform':platform.platform(),'input_sha256':hashlib.sha256(raw).hexdigest(),'actual_input_sha256':hashlib.sha256(actualraw).hexdigest(),
                    'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'wall_seconds':time.perf_counter()-start}}
    a.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k]for k in ['all_passed','menu_upper_count','joint_survivor_strict_upper','simple_strict_upper','gap_to_simple_upper']}))


if __name__=='__main__':main()
