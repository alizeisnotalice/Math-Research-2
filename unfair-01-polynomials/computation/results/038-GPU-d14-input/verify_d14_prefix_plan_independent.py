#!/usr/bin/env python3
"""Independent standard-library audit for d=14 prefix data.

The finite oracle is enumerated once.  The same production `validate` routine
then checks the accepted plan and every semantic mutation.  No generator or
worker code is imported and no MPS device is initialized.
"""
from __future__ import annotations
import argparse, copy, hashlib, itertools, json, os, sys, time
from collections import Counter
from pathlib import Path

OPTIONS={'0':(0,0),'1':(1,0),'2':(0,1)}
def digest(v): return hashlib.sha256(json.dumps(v,separators=(',',':')).encode()).hexdigest()
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def decode(code):
    a,b=[1],[1]
    for symbol in code:
        u,v=OPTIONS[symbol];a.append(u);b.append(v)
    return a,b
def coeff(a,b,k): return sum(a[i]*b[k-i] for i in range(max(0,k-len(b)+1),min(k,len(a)-1)+1))
def allowed(code,d):
    a,b=decode(code)
    return all(coeff(a,b,k) in (0,1) for k in range(d)) and sum(a[i]*b[d-i] for i in range(1,d))==0
def reasons(code,d):
    a,b=decode(code);product=[coeff(a,b,k) for k in range(2*d-1)];out=[]
    if max(product)>1: out.append('prefix_product_coefficient_gt_1')
    if any((product[d+i] if d+i<len(product) else 0)>=1 and a[i]+b[i]>0 for i in range(1,d)):
        out.append('forced_positive_fractional_cross_at_d_plus_i')
    return out
def oracle():
    raw=[''.join(chars) for chars in itertools.product('012',repeat=13)]
    valid=[code for code in raw if allowed(code,14)]
    survivors=[code for code in valid if not reasons(code,14)]
    starts=32768; assignments=[survivors[i%len(survivors)] for i in range(starts)];counts=Counter(assignments)
    reason_counts={name:sum(name in reasons(code,14) for code in valid) for name in ('prefix_product_coefficient_gt_1','forced_positive_fractional_cross_at_d_plus_i')}
    return {'raw':raw,'valid':valid,'survivors':survivors,'assignments':assignments,'counts':counts,'reasons':reason_counts}
def validate(plan, worker_path, prefix_path, worker_audit, O):
    """The sole production acceptance path for canonical data and mutations."""
    assert plan['schema']=='d14-necessary-prefix-plan-v1' and plan['status']=='cpu_generated_not_gpu_started'
    assert plan['parameters']=={'coverage_plan_starts':32768,'d':14,'m':60,'margin':.001,'minimum_starts':8192,'n':100,'seed':2026111615}
    assert plan['raw_ternary_count']==len(O['raw'])==3**13
    assert plan['valid_low_prefix_count']==len(O['valid']) and plan['survivor_count']==len(O['survivors'])
    assert plan['valid_low_prefix_codes_sha256']==digest(O['valid']) and plan['survivor_codes_sha256']==digest(O['survivors'])
    assert plan['survivor_codes']==O['survivors'] and plan['coverage_assignment_codes']==O['assignments']
    assert plan['per_prefix_assignment_counts']==[{'code':x,'count':O['counts'][x]} for x in O['survivors']]
    assert plan['minimum_starts_covers_pool'] is False and plan['coverage_plan_covers_every_survivor'] is True and plan['coverage_plan_is_8192_multiple'] is True
    assert min(O['counts'].values())==1 and max(O['counts'].values())==2 and sum(O['counts'].values())==32768
    worker_text=worker_path.read_text(); assert '1..14' in worker_text and 'd <= 14' in worker_text
    assert worker_audit['worker_sha256']==sha(worker_path)
    assert worker_audit['prefix_module_sha256']==sha(prefix_path)
    frozen=Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/ds_optimize/round033-gpu-d13/worker/gpu_factor_search_prefix.py')
    assert worker_audit['frozen_source_path']==str(frozen) and worker_audit['frozen_source_sha256']==sha(frozen)
    expected={'d':14,'raw_ternary_code_count':len(O['raw']),'valid_low_prefix_count':len(O['valid']),'survivor_count':len(O['survivors']),'valid_low_prefix_codes_sha256':digest(O['valid']),'survivor_codes_sha256':digest(O['survivors']),'strict_filter_rejected_by_reason_counts':O['reasons']}
    assert {key:worker_audit[key] for key in expected}==expected
    assert worker_audit['all_survivors_cpu_audited'] is True and worker_audit['float_arithmetic_used'] is False
    assert worker_audit['starts']==32768 and worker_audit['coverage_by_start_rows'] is True
def rejected(label, thunk):
    try: thunk()
    except (AssertionError, KeyError, TypeError): return label
    raise AssertionError('production validator accepted mutation: '+label)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--worker',type=Path,required=True);ap.add_argument('--prefix-module',type=Path,required=True);ap.add_argument('--worker-cpu-audit',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);x=ap.parse_args()
    assert os.path.ismount('/Volumes/data');begin=time.monotonic();plan=json.loads(x.input.read_text());wa=json.loads(x.worker_cpu_audit.read_text());O=oracle();validate(plan,x.worker,x.prefix_module,wa,O)
    p1=copy.deepcopy(plan);p1['survivor_codes'].pop()
    p2=copy.deepcopy(plan);p2['coverage_assignment_codes'][1]=O['survivors'][2]
    p3=copy.deepcopy(plan);p3['parameters']['coverage_plan_starts']=8192
    w1=copy.deepcopy(wa);w1['survivor_count']+=1
    w2=copy.deepcopy(wa);w2['worker_sha256']='0'*64
    w3=copy.deepcopy(wa);w3['prefix_module_sha256']='0'*64
    negatives=[
        rejected('plan_drop_survivor',lambda:validate(p1,x.worker,x.prefix_module,wa,O)),
        rejected('plan_wrong_rotation',lambda:validate(p2,x.worker,x.prefix_module,wa,O)),
        rejected('plan_false_coverage_bound',lambda:validate(p3,x.worker,x.prefix_module,wa,O)),
        rejected('worker_audit_count_mismatch',lambda:validate(plan,x.worker,x.prefix_module,w1,O)),
        rejected('worker_audit_worker_hash_mismatch',lambda:validate(plan,x.worker,x.prefix_module,w2,O)),
        rejected('worker_audit_prefix_hash_mismatch',lambda:validate(plan,x.worker,x.prefix_module,w3,O)),
    ]
    out={'status':'passed_independent_integer_d14_plan_audit_v2','d':14,'raw_ternary_count':len(O['raw']),'valid_low_prefix_count':len(O['valid']),'survivor_count':len(O['survivors']),'valid_low_prefix_codes_sha256':digest(O['valid']),'survivor_codes_sha256':digest(O['survivors']),'strict_filter_rejected_by_reason_counts':O['reasons'],'coverage_plan_starts':32768,'per_prefix_assignment_count_min':min(O['counts'].values()),'per_prefix_assignment_count_max':max(O['counts'].values()),'all_survivors_independently_checked':True,'all_coverage_rows_independently_checked':True,'production_validator_used_for_positive_and_negative_cases':True,'negative_tests_rejected':negatives,'source_sha256':{'plan':sha(x.input),'worker':sha(x.worker),'prefix_module':sha(x.prefix_module),'worker_cpu_audit':sha(x.worker_cpu_audit),'verifier':sha(Path(__file__))},'interpreter':sys.version,'elapsed_seconds':time.monotonic()-begin,'claim_boundary':'Independent finite integer audit only; no MPS/GPU was initialized and no factorization or exclusion is claimed.'}
    x.out.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps(out,sort_keys=True))
if __name__=='__main__': main()
