#!/usr/bin/env python3
"""Verify saved source/hash/strict-improvement chain; do not rerun geometry."""
import argparse
from fractions import Fraction as F
import hashlib
import json
import os
from pathlib import Path


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--pilot',type=Path,required=True)
    ap.add_argument('--continuation',type=Path,required=True);ap.add_argument('--seed',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True);args=ap.parse_args();assert os.path.ismount('/Volumes/data')
    pilot=json.loads((args.pilot/'result.json').read_text());spec=json.loads((args.pilot/'input.json').read_text())
    receipt=json.loads((args.pilot/'receipt.json').read_text())
    summary=json.loads((args.continuation/'run/continuation-summary.json').read_text())
    seed=json.loads((args.seed/'seed-deletion.json').read_text())
    module=args.pilot/'mass_polytope_variation.py';driver=args.continuation/'continue_certified_line.py'
    assert sha(module)==pilot['source_sha256']==summary['module_sha256']=='4ef61b815b9c4a9a8972cdaf31716ae7917b59243a29fb8697c8042137788246'
    assert sha(driver)==summary['driver_sha256']=='5e2f9f7dd0a3c97402251c751eb13307af25fd9182a1bdbbda67236bc78c8ed6'
    assert receipt['exit_status']==0 and receipt['owned_pid_absent']
    assert sha(args.pilot/'input.json')==pilot['input_sha256']==summary['original_spec_sha256']
    assert F(seed['deletions']['max_mass_min_lower'])<=F(pilot['seed_objective']['lower'])
    assert F(pilot['seed_objective']['upper'])<=F(seed['deletions']['max_mass_min_upper'])
    assert summary['status']=='FINITE_CALL_BUDGET_EXHAUSTED' and len(summary['steps'])==8
    assert summary['initial_result_sha256']==sha(args.pilot/'result.json')
    records=[];best=pilot['best_candidate'];previous=best['objective'];parent=sha(args.pilot/'result.json')
    line_candidates=len(pilot['candidates']);feasible=sum(r['exact_check']['feasible'] for r in pilot['candidates'])
    improvements=sum(bool(r.get('certified_strict_improvement')) for r in pilot['candidates'])
    pilot_record={'phase':'pilot','LP_calls':1,'line_points':len(pilot['candidates']),
                  'seed_objective':pilot['seed_objective'],'accepted_objective':previous,
                  'accepted_step':best['step'],'active_full_boxes':previous['full_active_boxes']}
    assert best['exact_check']['feasible'] and F(previous['lower'])>F(pilot['seed_objective']['upper'])
    for k,step in enumerate(summary['steps'],1):
        ip=args.continuation/f'run/input-{k:02}.json';op=args.continuation/f'run/result-{k:02}.json'
        inp=json.loads(ip.read_text());result=json.loads(op.read_text())
        assert sha(ip)==step['input_sha256']==result['input_sha256']
        assert sha(op)==step['output_sha256'] and result['source_sha256']==sha(module)
        assert parent==step['parent_result_sha256']==inp['continuation_parent']['result_sha256']
        assert inp['initial_masses']==best['masses']
        assert result['seed_objective']==previous
        assert step['accepted'] and step['exit_status']==0 and not step['stderr']
        best=result['best_candidate'];new=best['objective'];state=step['accepted_state']
        assert best['exact_check']['feasible'] and state['exact_check']['feasible']
        assert state['exact_check']['maximum_eligible_masks']==[spec['S_mask']]
        assert state['exact_check']['maximum_eligible_mass']=='1'
        assert len(state['connected_components'])==1 and state['positive_source_count']==21
        assert F(state['actual_nonterminal_margin'])>0 and F(state['minimum_mass'])>0
        assert F(new['lower'])>F(previous['upper'])
        assert F(step['strict_gain_lower'])==F(new['lower'])-F(previous['upper'])
        assert new==step['accepted_objective']
        line_candidates+=len(result['candidates']);feasible+=sum(r['exact_check']['feasible'] for r in result['candidates'])
        improvements+=sum(bool(r.get('certified_strict_improvement')) for r in result['candidates'])
        records.append({'call':k,'input_sha256':sha(ip),'parent_result_sha256':parent,'output_sha256':sha(op),
                        'line_points':len(result['candidates']),'accepted_step':best['step'],
                        'objective':new,'gain_lower':step['strict_gain_lower'],'state':state})
        previous=new;parent=sha(op)
    assert previous==summary['final_objective'] and best['masses']==summary['final_masses']
    payload={'round':64,'status':'SAVED HASH AND STRICT CERTIFICATE CHAIN VERIFIED; INDEPENDENT GEOMETRY AUDIT SEPARATE',
             'summary_source_sha256':sha(Path(__file__)),'pilot':pilot_record,'continuation_steps':records,
             'total_LP_calls':9,'total_line_candidates':line_candidates,'exactly_feasible_line_candidates':feasible,
             'strictly_improving_line_candidates':improvements,'accepted_updates_including_pilot':9,
             'termination':'FINITE_CALL_BUDGET_EXHAUSTED','converged':False,'local_optimality_claimed':False,
             'final_objective':summary['final_objective'],'final_masses':summary['final_masses'],
             'last_overlap_competitor_mass':records[-1]['state']['competitor_overlap_max_mass'],
             'root_seed_bits':seed['bits'],'root_seed_deletion_bounds':[seed['deletions']['max_mass_min_lower'],seed['deletions']['max_mass_min_upper']],
             'scope':'Same21positions, finite legal mass variation. Saved ledger/hash/interval comparisons only; no new geometry run by summarizer; no general C2 bound.'}
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(payload,indent=2)+'\n')
    print(json.dumps({k:payload[k] for k in ('total_LP_calls','total_line_candidates','exactly_feasible_line_candidates','strictly_improving_line_candidates','termination')}))


if __name__=='__main__':main()
