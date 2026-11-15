#!/usr/bin/env python3
"""Frozen geometry control for core's analytic five-source thin counterexample."""
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True)
    ap.add_argument('--budget',type=Path,required=True);ap.add_argument('--erosion',type=Path,required=True)
    ap.add_argument('--module',type=Path,required=True);ap.add_argument('--out-dir',type=Path,required=True)
    args=ap.parse_args();assert os.path.ismount('/Volumes/data')
    assert sha(args.input)=='b9e9888e8d9f2f12525b41dc85cfeda07f5554527fae838dc59ab4706fdd1940'
    assert sha(args.budget)=='211d33e21e47e6e14290e4ddbc3b3a36e247131de1bc53114091b7766d0fe8b6'
    assert sha(args.erosion)=='4dc4fb494736c40334359f9950f052ee2e4a53b503883ffbc3f9bc0da69b1361'
    args.out_dir.mkdir(parents=True,exist_ok=True);d=json.loads(args.input.read_text())
    adapter=dict(d);adapter['id']='r066-five-source-thin-shape-loss-control';adapter['original_exact_input_sha256']=sha(args.input)
    ip=args.out_dir/'input-adapter.json';ip.write_text(json.dumps(adapter,indent=2)+'\n')
    bp=args.out_dir/'budget.json';ep=args.out_dir/'erosion.json';receipts=[]
    commands=[[sys.executable,'-B',str(args.budget),'--module',str(args.module),'--r64-input',str(ip),'--out',str(bp)],
              [sys.executable,'-B',str(args.erosion),'--module',str(args.module),'--budget-result',str(bp),'--out',str(ep)]]
    for cmd in commands:
        t=time.monotonic();p=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
        stdout,stderr=p.communicate(timeout=60);ps=subprocess.run(['ps','-p',str(p.pid),'-o','pid='],capture_output=True,text=True)
        receipts.append({'command':cmd,'pid':p.pid,'exit_status':p.returncode,'wall_seconds':time.monotonic()-t,
                         'stdout':stdout,'stderr':stderr,'owned_pid_absent':ps.returncode!=0 and not ps.stdout.strip()})
        assert p.returncode==0 and receipts[-1]['owned_pid_absent']
    b=json.loads(bp.read_text())['cases'][0]['blocks'][0]
    D=1<<96;k=math.isqrt(3*D*D);rlo=F(k,D);rhi=F(k+1,D)
    quantities={'C_area':b['container_area'],'low_old_credit':b['low_old_credit'],
                'high_unpaid':b['high_unpaid_union'],'compensated_high':b['high_unpaid_minus_low_old_credit_normalized'],
                'shape_free_margin':b['universal_compensation_threshold'],'actual_shape_budget':b['analytic_envelope_budget']}
    checks={}
    for name,enclosure in quantities.items():
        a,c=map(F,d['quantities_as_rational_plus_coefficient_sqrt3'][name])
        lo,hi=(a+c*rlo,a+c*rhi) if c>=0 else (a+c*rhi,a+c*rlo)
        assert F(enclosure['lower'])<=hi and lo<=F(enclosure['upper']),name
        checks[name]={'closed_form_lower':str(lo),'closed_form_upper':str(hi),'geometry_interval_overlap':True}
    score=b['high_unpaid_minus_low_old_credit_normalized'];threshold=b['universal_compensation_threshold']
    assert F(score['lower'])>F(threshold['upper'])
    assert len(b['components'])==1 and not b['components'][0]['terminal_35_33']
    assert b['S_width_squared_over_L']==['9/16','0'] and b['budget_certifies_7over2']
    output={'status':'PASS','scope':'Five-source connected nonterminal thin example; strict failure of shape-free strong condition only, not fat-branch or general3.5 counterexample.',
            'python':sys.version,'runner_sha256':sha(Path(__file__)),'original_input_sha256':sha(args.input),
            'budget_sha256':sha(args.budget),'erosion_sha256':sha(args.erosion),'module_sha256':sha(args.module),
            'budget_output_sha256':sha(bp),'erosion_output_sha256':sha(ep),'receipts':receipts,
            'closed_form_checks':checks,'strict_strong_failure':True,'fat_branch':False}
    (args.out_dir/'receipt.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({'status':'PASS','strict_strong_failure':True,'fat_branch':False}))


if __name__=='__main__':main()
