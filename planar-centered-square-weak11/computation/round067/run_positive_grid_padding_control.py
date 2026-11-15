#!/usr/bin/env python3
"""Exactly one positive 9+1 grid-padding control; frozen R66 kernels unchanged."""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import time


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def load(name,p):
    spec=importlib.util.spec_from_file_location(name,p);mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod);return mod


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--budget',type=Path,required=True)
    ap.add_argument('--module',type=Path,required=True);ap.add_argument('--out-dir',type=Path,required=True)
    args=ap.parse_args();assert os.path.ismount('/Volumes/data')
    assert sha(args.budget)=='211d33e21e47e6e14290e4ddbc3b3a36e247131de1bc53114091b7766d0fe8b6'
    assert sha(args.module)=='4ef61b815b9c4a9a8972cdaf31716ae7917b59243a29fb8697c8042137788246'
    started=time.monotonic();mod=load('frozen_r064_geometry',args.module)
    budget=load('frozen_r066_budget',args.budget)
    d=F(1023,1024);eta=F(1,1<<20);epsilon=F(1,1024);theta=F(3,4)
    axis=[F(0),d/2,d];points=[(x,y) for y in axis for x in axis]+[(d+F(1,8),d)]
    masses=[eta]*10;masses[0]=F(1,5);masses[8]=F(4,5)-7*eta;masses[9]=epsilon
    S=(1<<9)-1;assert sum(masses[:9])==1 and min(masses)>0
    assert set(points[:9])=={(x,y) for x in axis for y in axis}
    # A separate 2^10 subset-to-closure enumeration checks the complete bbox
    # catalog against the frozen coordinate-interval method.
    catalog={}
    for subset in range(1,1<<10):
        ids=[i for i in range(10) if subset>>i&1]
        a=min(points[i][0] for i in ids);b=max(points[i][0] for i in ids)
        c=min(points[i][1] for i in ids);e=max(points[i][1] for i in ids)
        closure=sum(1<<i for i,(x,y) in enumerate(points) if a<=x<=b and c<=y<=e)
        mass=sum(masses[i] for i in range(10) if closure>>i&1);span2=max(b-a,e-c)**2
        catalog[closure]={'mask':closure,'bbox':list(map(str,(a,b,c,e))),
                          'mass':str(mass),'span_squared':str(span2),'strict_eligible':span2<mass}
    frozen_boxes=mod.closures(points);assert set(catalog)=={b['mask'] for b in frozen_boxes}
    eligible=[r for r in catalog.values() if r['strict_eligible']]
    L=max(F(r['mass']) for r in eligible);maxima=sorted(r['mask'] for r in eligible if F(r['mass'])==L)
    assert L==1 and maxima==[S]
    connected=budget.components(points,masses,[r['mask'] for r in eligible])
    assert len(connected)==1 and not connected[0]['terminal_35_33']
    assert 16*d*d>9 and (d/2)**2<theta
    # Exact product-grid overlap plus direct integer union equalities certify
    # that the uniform source-neighborhood union is exactly old C_theta.
    bits=72;D=1<<bits;ri,ro=budget.root_radius(theta,bits)
    uniform_checks=[]
    for radius in (ri,ro):
        squares=[]
        for i,(x,y) in enumerate(points[:9]):
            xx=x*D;yy=y*D;assert xx.denominator==yy.denominator==1
            squares.append((int(xx)-radius,int(xx)+radius,int(yy)-radius,int(yy)+radius,1<<i))
        endpoint=d*D;assert endpoint.denominator==1
        C=(-radius,int(endpoint)+radius,-radius,int(endpoint)+radius,S)
        outside=mod.difference_area(squares,[C]);missing=mod.difference_area([C],squares)
        assert outside==missing==0
        uniform_checks.append({'radius_numerator':str(radius),'union_minus_container_area_numerator':str(outside),
                               'container_minus_union_area_numerator':str(missing)})
    data={'schema':'r067-positive-grid-padding-control-v1','id':'r067-positive-nine-plus-one-uniform-neighborhood-obstruction',
          'parameters':{'d':str(d),'eta':str(eta),'epsilon':str(epsilon),'theta':str(theta)},
          'points':[[str(x),str(y)] for x,y in points],'initial_masses':list(map(str,masses)),
          'S_mask':S,'S_indices':list(range(9)),'source_roles':['S']*9+['outside'],
          'theta':str(theta),'S_mass':'1','enumerated_nonempty_source_subsets':1023,
          'closed_bbox_catalog':[catalog[s] for s in sorted(catalog)],'closed_bbox_count':len(catalog),
          'strict_eligible_count':len(eligible),'maximum_eligible_mass':str(L),'all_maximum_masks':maxima,
          'positive_source_count':10,'components':connected,'S_width_squared_over_L':[str(d*d)]*2,
          'uniform_support_neighborhood_union_equals_C':True,
          'grid_gap_squared_less_than_theta_L':str((d/2)**2-theta),
          'uniform_neighborhood_directed_union_checks':uniform_checks,
          'builder_sha256':sha(Path(__file__)),
          'scope':'Single prescribed positive grid-padding input. Tests uniform-radius source-support union only; no LP/random extension or new mass-sensitive state claim.'}
    args.out_dir.mkdir(parents=True,exist_ok=True);ip=args.out_dir/'input.json'
    ip.write_text(json.dumps(data,indent=2)+'\n');op=args.out_dir/'budget.json'
    command=[sys.executable,'-B',str(args.budget),'--module',str(args.module),'--r64-input',str(ip),'--out',str(op)]
    t=time.monotonic();proc=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    try:stdout,stderr=proc.communicate(timeout=60)
    except subprocess.TimeoutExpired:
        proc.kill();proc.communicate();raise
    ps=subprocess.run(['ps','-p',str(proc.pid),'-o','pid='],capture_output=True,text=True)
    execution={'command':command,'pid':proc.pid,'exit_status':proc.returncode,'wall_seconds':time.monotonic()-t,
               'stdout':stdout,'stderr':stderr,'owned_pid_absent':ps.returncode!=0 and not ps.stdout.strip()}
    assert proc.returncode==0 and execution['owned_pid_absent']
    result=json.loads(op.read_text());assert result['status']=='PASS'
    b=result['cases'][0]['blocks'][0];assert b['S_mask']==S and b['closed_bbox_count']==len(catalog)
    assert b['strict_eligible_count']==len(eligible)
    assert F(b['analytic_envelope_budget']['lower'])>F(7,2)
    assert F(b['true_deletion']['upper'])<F(7,2)
    receipt={'status':'PASS_POSITIVE_GRID_PADDING_FIXED_THETA_ENVELOPE_OBSTRUCTION','python':sys.version,
             'builder_sha256':sha(Path(__file__)),'budget_source_sha256':sha(args.budget),'geometry_source_sha256':sha(args.module),
             'input_sha256':sha(ip),'output_sha256':sha(op),'execution':execution,'total_wall_seconds':time.monotonic()-started,
             'closed_bbox_count':len(catalog),'strict_eligible_count':len(eligible),'unique_maximum_mask':S,
             'uniform_support_neighborhood_union_equals_C':True,'production_budget':b['analytic_envelope_budget'],
             'production_actual_deletion':b['true_deletion'],'scope':data['scope']}
    (args.out_dir/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':receipt['status'],'closed_boxes':len(catalog),'eligible':len(eligible),
                      'budget_lower_display':float(F(b['analytic_envelope_budget']['lower'])),
                      'deletion_upper_display':float(F(b['true_deletion']['upper']))}))


if __name__=='__main__':main()
