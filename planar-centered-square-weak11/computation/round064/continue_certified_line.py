#!/usr/bin/env python3
"""Bounded same-position continuation; accept only successive strict certificates."""
import argparse
import copy
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


def write(path,data):
    path.parent.mkdir(parents=True,exist_ok=True);tmp=path.with_suffix(path.suffix+'.partial')
    tmp.write_text(json.dumps(data,indent=2)+'\n');tmp.replace(path)


def state_record(mod,spec,masses):
    p=[tuple(map(F,q)) for q in spec['points']];boxes=mod.closures(p)
    check=mod.exact_check(masses,spec,mod.constraints(spec,boxes),boxes)
    assert check['feasible'],check
    eligible=[];components=[1<<i for i in range(len(p))]
    for b in boxes:
        mass=sum(masses[i] for i in b['ids'])
        if mass<=b['span2']:continue
        mask=b['mask'];joined=mask
        touched=[c for c in components if c&mask]
        for c in touched:joined |= c
        components=[c for c in components if not c&mask]+[joined]
        eligible.append({'mask':mask,'rank':mask.bit_count(),'mass':str(mass),
                         'span_squared':str(b['span2'])})
    d2=max(max(q[a] for q in p)-min(q[a] for q in p) for a in (0,1))**2
    margin=d2-F(35,33)*sum(masses)
    assert len(components)==1 and margin>0
    floor=F(spec['mass_constraints']['positive_floor_for_initial_chamber'])
    competitors=[F(r['mass']) for r in eligible if r['mask']!=spec['S_mask']]
    overlap_competitors=[F(r['mass']) for r in eligible
                         if r['mask']!=spec['S_mask'] and r['mask']&spec['S_mask']]
    return {'exact_check':check,'eligible_closures':eligible,'connected_components':components,
            'positive_source_count':sum(m>0 for m in masses),'actual_nonterminal_margin':str(margin),
            'minimum_mass':str(min(masses)),'floor_active_indices':[i for i,m in enumerate(masses) if m==floor],
            'eligible_competitor_max_mass_excluding_exact_S':str(max(competitors,default=F(0))),
            'competitor_overlap_max_mass':str(max(overlap_competitors,default=F(0))),
            'competitor_scope':'Maximum strict eligible bbox-closure mass except exact S, not an independently reconstructed positive-area N-state catalog.'}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True)
    ap.add_argument('--initial-result',type=Path,required=True);ap.add_argument('--module',type=Path,required=True)
    ap.add_argument('--out-dir',type=Path,required=True);ap.add_argument('--max-calls',type=int,default=8)
    ap.add_argument('--seconds-per-call',type=int,default=60);ap.add_argument('--total-seconds',type=int,default=600)
    args=ap.parse_args();assert os.path.ismount('/Volumes/data')
    assert 1<=args.max_calls<=8 and args.seconds_per_call>0 and args.total_seconds>0
    expected='4ef61b815b9c4a9a8972cdaf31716ae7917b59243a29fb8697c8042137788246'
    assert sha(args.module)==expected,'Frozen LP module identity changed'
    spec=json.loads(args.input.read_text());initial=json.loads(args.initial_result.read_text())
    assert initial['input_sha256']==sha(args.input) and initial['source_sha256']==expected
    best=initial['best_candidate'];assert best['exact_check']['feasible']
    assert F(best['objective']['lower'])>F(initial['seed_objective']['upper'])
    masses=list(map(F,best['masses']));objective=best['objective']
    mspec=importlib.util.spec_from_file_location('frozen_lp_module',args.module)
    mod=importlib.util.module_from_spec(mspec);mspec.loader.exec_module(mod)
    args.out_dir.mkdir(parents=True,exist_ok=True);summary_path=args.out_dir/'continuation-summary.json'
    summary={'schema':'r064-bounded-same-position-continuation-v1','status':'RUNNING',
             'driver_sha256':sha(Path(__file__)),'module_sha256':expected,
             'original_spec_sha256':sha(args.input),'initial_result_sha256':sha(args.initial_result),
             'original_seed_objective':initial['seed_objective'],'starting_objective':objective,
             'starting_state':state_record(mod,spec,masses),'steps':[],
             'budgets':{'max_calls':args.max_calls,'seconds_per_call':args.seconds_per_call,'total_seconds':args.total_seconds},
             'scope':'Successive strict full-halo certificates on same fixed positions; no local optimum or class upper bound claim.'}
    started=time.monotonic();trust=F(1,16);nogain=0;parent_sha=sha(args.initial_result)
    write(summary_path,summary)
    for call in range(args.max_calls):
        remaining=args.total_seconds-(time.monotonic()-started)
        if remaining<=0:summary['status']='TOTAL_TIME_BUDGET_REACHED';break
        ss=copy.deepcopy(spec);ss['id']=spec['id']+f'-continuation-{call+1:02}'
        ss['initial_masses']=list(map(str,masses))
        if 'actual_seed_checks' in ss:ss['historical_original_seed_checks']=ss.pop('actual_seed_checks')
        ss['status']='Exact previous certified winner used as next initial masses; current checks in result.'
        ss['continuation_parent']={'result_sha256':parent_sha,'original_input_sha256':sha(args.input)}
        ip=args.out_dir/f'input-{call+1:02}.json';op=args.out_dir/f'result-{call+1:02}.json'
        write(ip,ss)
        command=[sys.executable,'-B',str(args.module),'--input',str(ip),'--out',str(op),'--trust-radius',str(trust)]
        record={'call':call+1,'command':command,'trust_radius':str(trust),'input_sha256':sha(ip),'parent_result_sha256':parent_sha}
        t=time.monotonic()
        try:
            proc=subprocess.run(command,capture_output=True,text=True,timeout=min(args.seconds_per_call,remaining))
            record.update({'exit_status':proc.returncode,'wall_seconds':time.monotonic()-t,'stdout':proc.stdout,'stderr':proc.stderr})
            assert proc.returncode==0,'LP child failed'
            result=json.loads(op.read_text());record['output_sha256']=sha(op)
            assert result['input_sha256']==sha(ip) and result['source_sha256']==expected
            assert F(result['seed_objective']['lower'])==F(objective['lower'])
            assert F(result['seed_objective']['upper'])==F(objective['upper'])
            candidate=result.get('best_candidate')
            accepted=bool(candidate and candidate['exact_check']['feasible'] and
                          F(candidate['objective']['lower'])>F(objective['upper']))
            record['accepted']=accepted
            if accepted:
                assert candidate['certified_strict_improvement']
                nm=list(map(F,candidate['masses']));record['accepted_state']=state_record(mod,spec,nm)
                record['strict_gain_lower']=str(F(candidate['objective']['lower'])-F(objective['upper']))
                masses=nm;objective=candidate['objective'];parent_sha=sha(op)
                record['accepted_objective']=objective;trust=F(1,16);nogain=0
            else:trust/=2;nogain+=1
            summary['steps'].append(record)
            summary['final_masses']=list(map(str,masses));summary['final_objective']=objective
            summary['wall_seconds']=time.monotonic()-started;write(summary_path,summary)
            print(json.dumps({'call':call+1,'accepted':accepted,'objective_lower':objective['lower'],'no_gain_streak':nogain}),flush=True)
            if nogain>=2:summary['status']='TWO_NO_GAIN_CALLS_STOPPED';break
        except Exception as exc:
            record['error']=repr(exc);summary['steps'].append(record);summary['status']='FAILED_RETAIN_PREVIOUS_CERTIFIED_WINNER'
            write(summary_path,summary);raise
    else:summary['status']='FINITE_CALL_BUDGET_EXHAUSTED'
    summary['wall_seconds']=time.monotonic()-started;write(summary_path,summary)


if __name__=='__main__':main()
