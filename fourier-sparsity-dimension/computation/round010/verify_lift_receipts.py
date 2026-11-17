#!/usr/bin/env python3
"""Bind all lift samples to sources, independently check their parameters,
and compare every field of the replicated 320-sample local/remote job."""
import argparse,hashlib,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args();base=Path(__file__).parent
read=lambda name:json.loads((base/name).read_text())
remote=read('remote-lift-execution.json');local=read('local-lift-execution.json');seed=read('independent-seed-verification.json')
assert remote['ssh_exit_code']==local['exit_code']==0
r=remote['result'];assert r['status']=='passed'
source_hash=hashlib.sha256((base/'cddl_fact22_lift_verify.cpp').read_bytes()).hexdigest()
assert r['source_sha256']==local['source_sha256']==source_hash
truth=int(seed['truthmask_hex'],16);assert seed['status']=='passed'
seen=set();total=0;coefficients=0
for job in r['jobs']:
 assert job['exit_code']==0
 result=job['result'];assert result['status']=='passed';assert result['seed'] not in seen;seen.add(result['seed'])
 assert result['degree14_modulus_hex']==seed['modulus_hex'] and result['alpha_hex']==seed['alpha_hex']
 assert int(result['seed_truthmask_hex_le'],16)==truth
 assert result['seed_walsh_histogram']==seed['walsh_histogram'] and result['seed_anf_degree']==4
 assert len(result['samples'])==320
 counts={}
 for s in result['samples']:
  ell=s['l'];Q=1<<ell;counts[ell]=counts.get(ell,0)+1
  assert s['variables']==14+2*ell and s['q']==128*Q
  tables=s['triangular_prefix_truth_tables'];assert len(tables)==ell
  assert all(0<=t<(1<<(1<<j)) for j,t in enumerate(tables))
  pi=[w^sum(((tables[j]>>(w&((1<<j)-1)))&1)<<j for j in range(ell)) for w in range(Q)]
  assert sorted(pi)==list(range(Q))
  assert len(s['lin_by_w'])==len(s['constant_by_w'])==Q
  assert all(0<=x<1<<14 for x in s['lin_by_w']) and all(x in (0,1) for x in s['constant_by_w'])
  assert all(s[key] is True for key in ('pi_bijective','all_walsh_abs_q','direct_walsh_matches_formula','one_bit_flip_breaks_flatness'))
  total+=1;coefficients+=1<<s['variables']
 assert counts=={i:64 for i in range(5)}
assert len(r['jobs'])==14
left=dict(local['result']);right=dict(r['jobs'][0]['result']);left.pop('elapsed_seconds');right.pop('elapsed_seconds');assert left==right
output={'status':'passed','arithmetic':'exact integer/F2','remote_jobs':14,'remote_samples':total,
 'direct_walsh_coefficients_checked_remotely':coefficients,'cross_machine_samples_compared':len(left['samples']),
 'independent_seed_truth_table_matches':True,'source_sha256_bound':source_hash,
 'peak_remote_task_cpu_percent':max(s['percent_sum'] for s in r['cpu_samples']),
 'remote_whole_machine_top':r['whole_machine_top'],'remote_worker_elapsed_seconds':[j['result']['elapsed_seconds'] for j in r['jobs']],
 'scope':'finite Walsh/formula checks, not a new nonweak-normality enumeration; all-parameter legality follows separate proof',
 'inputs_sha256':{name:hashlib.sha256((base/name).read_bytes()).hexdigest() for name in ['remote-lift-execution.json','local-lift-execution.json','independent-seed-verification.json']},
 'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
a.out.write_text(json.dumps(output,indent=2)+'\n');print(json.dumps(output))
