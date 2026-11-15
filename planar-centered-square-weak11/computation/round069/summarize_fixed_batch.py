#!/usr/bin/env python3
"""Verify the saved eight-job source/receipt chain and summarize exact intervals."""
import argparse
from fractions import Fraction as F
import hashlib
import json
import os
from pathlib import Path
import sys


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--archive-root',type=Path,required=True)
    ap.add_argument('--spec',type=Path,required=True);ap.add_argument('--results',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True);args=ap.parse_args();assert os.path.ismount('/Volumes/data')
    sp=args.spec if args.spec.is_absolute() else args.archive_root/args.spec
    spec=json.loads(sp.read_text());receipt=json.loads((args.results/'batch-receipt.json').read_text())
    assert receipt['status']=='PASS' and len(receipt['jobs'])==len(spec['jobs'])==8
    rows=[];central_positive=[]
    for job in spec['jobs']:
        p=args.results/f"{job['job_id']}.json";d=json.loads(p.read_text())
        r=next(j for j in receipt['jobs'] if j['job_id']==job['job_id'])
        assert r['exit_status']==0 and r['owned_pid_absent'] and not r['stderr']
        assert sha(p)==r['output_sha256'] and sha(args.archive_root/job['input'])==r['input_sha256']==d['input_sha256']==job['input_sha256']
        assert d['source_sha256']==r['source_sha256']==receipt['source_sha256']
        assert d['spec_sha256']==r['spec_sha256']==sha(sp)
        assert d['status']=='PASS' and d['central_residual_halo_geometrically_certified_empty']
        regions=d['regions'];net=d['mixed_mass_gap_diagnostics']
        assert F(regions['R']['central']['upper'])==0
        for part in ['global','central','peripheral']:
            lo=F(net[part]['P_minus_credit']['normalized_lower']);hi=F(net[part]['P_minus_credit']['normalized_upper'])
            assert net[part]['net_at_most_L_over_5_certified']==(hi<=F(1,5))
            assert net[part]['net_strictly_above_L_over_5_certified']==(lo>F(1,5))
        if F(regions['P_pure_mixed']['central']['lower'])>0:central_positive.append(job['job_id'])
        rows.append({'id':d['id'],'source_input_id':d['source_input_id'],'S_mask':d['S_mask'],'L':d['L'],
          'all_maximum_masks':d['all_maximum_masks'],'source_count':d['source_count'],'S_rank':d['S_rank'],
          'W':regions['W'],'Wprime':regions['Wprime'],'Delta':regions['Delta'],'I':regions['I'],
          'credit':regions['internal_old_credit'],'P':regions['P_pure_mixed'],
          'net_mixed':net,'W_old_surplus':regions['W_old_surplus'],'Wprime_old_surplus':regions['Wprime_old_surplus'],
          'Wprime_surplus_outside_I':regions['Wprime_surplus_outside_I'],'K':d['K'],
          'R_K_upper_zero':True,'input_sha256':d['input_sha256'],'output_sha256':sha(p),
          'domain':{'component_masks':d['eligible_component_masks'],'nonterminal_margin':d['whole_nonterminal_margin'],
                    'width_squared_over_L':d['S_width_squared_over_L']}})
    payload={'round':69,'status':'SAVED OUTPUT/SOURCE/PID RECEIPT CHAIN VERIFIED; NO HALO RERUN',
      'python':sys.version,'source_sha256':sha(Path(__file__)),'producer_source_sha256':receipt['source_sha256'],
      'spec_sha256':sha(sp),'geometry_jobs':8,'historical_source_inputs':5,'historical_maximum_states':6,'structural_controls':2,
      'geometry_child_exit0_and_PID_absent':True,'batch_wall_seconds':receipt['wall_seconds'],
      'sampled_max_rss_kib':max((s['rss_kib'] for r in receipt['jobs'] for s in r['samples']),default=None),
      'R_K_upper_zero_jobs':8,'central_pure_mixed_certified_positive_jobs':central_positive,
      'global_net_at_most_L_over_5_jobs':sum(r['net_mixed']['global']['net_at_most_L_over_5_certified'] for r in rows),
      'global_net_strictly_above_L_over_5_jobs':sum(r['net_mixed']['global']['net_strictly_above_L_over_5_certified'] for r in rows),
      'rows':rows,'scope':'Exactly specified historical and structural jobs. Passing the .2 diagnostic on these inputs is not a general claim. Central common intervals are compatibility checks, while residual exclusion is certified by upper0 here and proved generally in R68.'}
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(payload,indent=2)+'\n')
    print(json.dumps({'jobs':8,'net_pass':payload['global_net_at_most_L_over_5_jobs'],
                      'net_strict_fail':payload['global_net_strictly_above_L_over_5_jobs'],
                      'central_P_positive_jobs':central_positive}))


if __name__=='__main__':main()
