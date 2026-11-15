#!/usr/bin/env python3
"""One deterministic fat connected counterexample to fixed-theta ENVELOPE.

This is not a counterexample to the true deletion target. Frozen production
sources are reused unmodified; exact analytic lower bounds are also saved.
"""
import argparse
from fractions import Fraction as F
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--budget',type=Path,required=True)
    ap.add_argument('--erosion',type=Path,required=True);ap.add_argument('--module',type=Path,required=True)
    ap.add_argument('--out-dir',type=Path,required=True);args=ap.parse_args();assert os.path.ismount('/Volumes/data')
    assert sha(args.budget)=='211d33e21e47e6e14290e4ddbc3b3a36e247131de1bc53114091b7766d0fe8b6'
    assert sha(args.erosion)=='4dc4fb494736c40334359f9950f052ee2e4a53b503883ffbc3f9bc0da69b1361'
    assert sha(args.module)=='4ef61b815b9c4a9a8972cdaf31716ae7917b59243a29fb8697c8042137788246'
    args.out_dir.mkdir(parents=True,exist_ok=True)
    alpha=F(1023,1024);epsilon=F(1,1024);delta=F(1,8)
    points=[(F(0),F(0)),(alpha,alpha),(alpha+delta,alpha)]
    masses=[F(1,5),F(4,5),epsilon];subsets=[]
    for mask in range(1,8):
        ids=[i for i in range(3) if mask>>i&1];mass=sum(masses[i] for i in ids)
        span2=max(max(points[i][a] for i in ids)-min(points[i][a] for i in ids) for a in (0,1))**2
        subsets.append({'mask':mask,'mass':str(mass),'span_squared':str(span2),'strict_eligible':span2<mass})
    eligible=[r for r in subsets if r['strict_eligible']];maximum=max(F(r['mass']) for r in eligible)
    assert maximum==1 and [r['mask'] for r in eligible if F(r['mass'])==maximum]==[3]
    assert 16*alpha*alpha>9 and 33*(alpha+delta)**2>35*sum(masses)
    D=1<<128;k3=math.isqrt(3*D*D);k15=math.isqrt(3*D*D//20)
    assert k3*k3<=3*D*D<(k3+1)**2
    assert 20*k15*k15<=3*D*D<20*(k15+1)**2
    r3lo=F(k3,D);r3hi=F(k3+1,D);r15lo=F(k15,D);r15hi=F(k15+1,D)
    heavy_out_lo=F(33,80)-r15hi;heavy_out_hi=F(33,80)-r15lo
    C_lo=alpha*alpha+F(3,4)+alpha*r3lo;C_hi=alpha*alpha+F(3,4)+alpha*r3hi
    analytic_B_lo=C_lo-epsilon+heavy_out_lo
    assert analytic_B_lo>F(7,2)
    data={'schema':'r066-three-source-fat-envelope-counterexample-v1',
          'id':'r066-three-source-fat-fixedtheta-envelope-obstruction',
          'points':[[str(x),str(y)] for x,y in points],'initial_masses':list(map(str,masses)),
          'S_mask':3,'S_indices':[0,1],'theta':'3/4','all_nonempty_subsets':subsets,
          'maximum_eligible_mass':'1','unique_maximum_mask':3,'eligible_connected_by_masks':[3,6],
          'fat_width_squared_over_L':[str(alpha*alpha)]*2,
          'whole_nonterminal_margin':str((alpha+delta)**2-F(35,33)*sum(masses)),
          'analytic_heavy_singleton_outside_C':{'lower':str(heavy_out_lo),'upper':str(heavy_out_hi),'formula':'33/80-sqrt(3/20)'},
          'analytic_C_area':{'lower':str(C_lo),'upper':str(C_hi),'formula':'alpha^2+3/4+alpha sqrt3'},
          'analytic_B_lower':str(analytic_B_lo),'analytic_B_minus_7over2_lower':str(analytic_B_lo-F(7,2)),
          'scope':'Unique eligible maximum S is fat; input connected/nonterminal. Counterexample only to fixedtheta=.75 container-envelope/strong compensation selection, not actual deletion<=3.5.',
          'builder_sha256':sha(Path(__file__))}
    ip=args.out_dir/'input.json';ip.write_text(json.dumps(data,indent=2)+'\n')
    bp=args.out_dir/'budget.json';ep=args.out_dir/'erosion.json';receipts=[]
    commands=[[sys.executable,'-B',str(args.budget),'--module',str(args.module),'--r64-input',str(ip),'--out',str(bp)],
              [sys.executable,'-B',str(args.erosion),'--module',str(args.module),'--budget-result',str(bp),'--out',str(ep)]]
    for cmd in commands:
        t=time.monotonic();p=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
        try:stdout,stderr=p.communicate(timeout=60)
        except subprocess.TimeoutExpired:
            p.kill();p.communicate();raise
        ps=subprocess.run(['ps','-p',str(p.pid),'-o','pid='],capture_output=True,text=True)
        receipts.append({'command':cmd,'pid':p.pid,'exit_status':p.returncode,'wall_seconds':time.monotonic()-t,
                         'stdout':stdout,'stderr':stderr,'owned_pid_absent':ps.returncode!=0 and not ps.stdout.strip()})
        assert p.returncode==0 and receipts[-1]['owned_pid_absent']
    b=json.loads(bp.read_text())['cases'][0]['blocks'][0]
    score=b['high_unpaid_minus_low_old_credit_normalized'];threshold=b['universal_compensation_threshold']
    assert F(b['analytic_envelope_budget']['lower'])>F(7,2)
    assert F(score['lower'])>F(threshold['upper'])
    assert F(b['true_deletion']['upper'])<F(7,2)
    assert len(b['components'])==1 and not b['components'][0]['terminal_35_33']
    assert all(16*F(v)>9 for v in b['S_width_squared_over_L'])
    receipt={'status':'PASS_STRICT_FAT_FIXED_THETA_ENVELOPE_COUNTEREXAMPLE','python':sys.version,
             'builder_sha256':sha(Path(__file__)),'input_sha256':sha(ip),'budget_source_sha256':sha(args.budget),
             'erosion_source_sha256':sha(args.erosion),'module_sha256':sha(args.module),
             'budget_output_sha256':sha(bp),'erosion_output_sha256':sha(ep),'receipts':receipts,
             'analytic_lower_bound':str(analytic_B_lo),'production_budget':b['analytic_envelope_budget'],
             'production_actual_deletion':b['true_deletion'],'fat_branch':True,'connected':True,
             'scope':data['scope']}
    (args.out_dir/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':receipt['status'],'analytic_B_lower_display':float(analytic_B_lo),
                      'production_B_lower_display':float(F(b['analytic_envelope_budget']['lower'])),
                      'actual_deletion_upper_display':float(F(b['true_deletion']['upper']))}))


if __name__=='__main__':main()
