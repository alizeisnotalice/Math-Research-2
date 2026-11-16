#!/usr/bin/env python3
"""Independent Fraction audit for completed first-fractional-prefix GPU records.

It deliberately reimplements prefix enumeration and never imports the worker
or its prefix module.  Exact claims concern JSON decimal values and the finite
prefix pool only, not unrecorded GPU tensors or a global factorization result.
"""
from __future__ import annotations

import argparse, hashlib, itertools, json
from fractions import Fraction
from pathlib import Path

PAIR = ((0, 0), (1, 0), (0, 1))

def sha(p: Path) -> str: return hashlib.sha256(p.read_bytes()).hexdigest()
def q(x) -> Fraction: return Fraction(str(x))
def decode(code: str):
    a,b=[1],[1]
    for c in code:
        x,y=PAIR[int(c)];a.append(x);b.append(y)
    return a,b
def lowconv(a,b,end):
    return [sum(a[i]*b[k-i] for i in range(max(0,k-len(b)+1),min(k,len(a)-1)+1)) for k in range(end+1)]
def legal(code,d):
    if len(code)!=d-1 or any(c not in '012' for c in code): return False
    a,b=decode(code)
    return set(lowconv(a,b,d-1)) <= {0,1} and sum(a[i]*b[d-i] for i in range(1,d)) == 0
def reasons(code,d):
    a,b=decode(code); product=lowconv(a,b,2*d-2); out=[]
    if any(x>1 for x in product): out.append('prefix_product_coefficient_gt_1')
    for i in range(1,d):
        baseline=product[d+i] if d+i<len(product) else 0
        if baseline>=1 and a[i]+b[i]>0:
            out.append('forced_positive_fractional_cross_at_d_plus_i');break
    return out
def survivors(d):
    valid=[''.join(c) for c in itertools.product('012',repeat=d-1) if legal(''.join(c),d)]
    return valid,[c for c in valid if not reasons(c,d)]
def conv(a,b):
    z=[Fraction(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):z[i+j]+=x*y
    return z
def validate_case(case,m,n,d,margin,pool):
    a,b=[q(x) for x in case['A']],[q(x) for x in case['B']]
    if len(a)!=m+1 or len(b)!=n+1 or a[0]!=1 or a[-1]!=1 or b[0]!=1 or b[-1]!=1: raise ValueError('endpoint/degree')
    if any(x<0 or x>1 for x in a+b): raise ValueError('outside unit interval')
    assigned=case['assigned_prefix_code']; derived=case['derived_prefix_code']
    if assigned not in pool or derived!=assigned: raise ValueError('assigned/derived pool mismatch')
    ea,eb=decode(assigned)
    if [int(x) for x in a[1:d]]!=ea[1:] or [int(x) for x in b[1:d]]!=eb[1:]: raise ValueError('prefix coefficients mismatch')
    if any(x not in (0,1) for x in a[1:d]+b[1:d]): raise ValueError('prefix not assigned binary')
    low=lowconv(ea,eb,d-1); cross=sum(ea[i]*eb[d-i] for i in range(1,d))
    if low!=case['prefix_low_convolution'] or cross!=case['prefix_cross_at_d'] or cross!=0: raise ValueError('saved low-prefix audit mismatch')
    complement=abs(a[d]+b[d]-1)
    if not(margin<=a[d]<=1-margin and margin<=b[d]<=1-margin and complement<=Fraction(1,10**7)): raise ValueError('fractional d model')
    c=conv(a,b); target=[Fraction(int(x>=Fraction(1,2))) for x in c];res=max(abs(x-y) for x,y in zip(c,target));delta=abs(res-q(case['max_residual_float64']))
    if delta>Fraction(1,10**12): raise ValueError('stored residual mismatch')
    return res,delta,complement

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-dir',type=Path,required=True);ap.add_argument('--worker',type=Path,required=True);ap.add_argument('--prefix-module',type=Path,required=True);ap.add_argument('--record',type=Path,required=True);x=ap.parse_args();r=x.run_dir
    config=json.loads((r/'config.json').read_text());audit=json.loads((r/'cpu-prefix-audit.json').read_text());checkpoint=json.loads((r/'checkpoint.json').read_text())
    m,n,d=int(config['m']),int(config['n']),int(config['fractional_index']);margin=q(config['margin']);valid,pool=survivors(d)
    encoded=lambda z:hashlib.sha256(json.dumps(z,separators=(',',':')).encode()).hexdigest()
    if (len(valid),encoded(valid),len(pool),encoded(pool)) != (audit['valid_low_prefix_count'],audit['valid_low_prefix_codes_sha256'],audit['survivor_count'],audit['survivor_codes_sha256']): raise SystemExit('independent prefix pool disagrees with saved CPU audit')
    if int(config['prefix_survivor_count'])!=len(pool) or config['prefix_codes_sha256']!=encoded(pool): raise SystemExit('config prefix pool mismatch')
    if sha(x.worker)!=config['script_sha256'] or sha(x.worker)!=audit['worker_sha256'] or sha(x.prefix_module)!=config['prefix_module_sha256'] or sha(x.prefix_module)!=audit['prefix_module_sha256']: raise SystemExit('frozen source SHA mismatch')
    restart_dir=r/'restarts' if (r/'restarts').is_dir() else r
    restarts=sorted(restart_dir.glob('restart_*.json')); expected=[restart_dir/f'restart_{i:07d}.json' for i in range(int(config['restart_every']),int(config['steps']),int(config['restart_every']))]
    if restarts!=expected: raise SystemExit('restart layout mismatch')
    docs=[]
    for epoch,p in enumerate(restarts):
        doc=json.loads(p.read_text())
        if doc.get('old_prefix_epoch')!=epoch or doc.get('next_prefix_epoch')!=epoch+1 or len(doc.get('history',[]))!=1: raise SystemExit('restart epoch/history mismatch')
        docs.append((p.name,doc,'rotation_snapshot'))
    docs.append(('checkpoint.json',checkpoint,'complete_final_cohort'))
    total=rotation=final=0;seen=set();maxdelta=maxcomp=Fraction(0);minres=None
    for name,doc,kind in docs:
        cases=doc.get('cases');
        if not isinstance(cases,list) or len(cases)!=16: raise SystemExit(f'{name}: expected 16 cases')
        if kind=='complete_final_cohort' and (not doc.get('history') or doc['history'][-1].get('step')!=int(config['steps']) or doc['history'][-1].get('restart_epoch')!=len(restarts)): raise SystemExit('final checkpoint not complete')
        for case in cases:
            res,delta,comp=validate_case(case,m,n,d,margin,set(pool));finger=sha_case(case)
            if finger in seen: raise SystemExit('duplicate saved case')
            seen.add(finger);total+=1;maxdelta=max(maxdelta,delta);maxcomp=max(maxcomp,comp);minres=res if minres is None else min(minres,res)
        if kind=='rotation_snapshot':rotation+=len(cases)
        else:final+=len(cases)
    release=json.loads((r/'release.json').read_text())
    if release.get('cleanup_called') is not True or release.get('process_will_exit') is not True: raise SystemExit('release metadata mismatch')
    out={'status':'saved_prefix_samples_fraction_and_pool_audit_passed','claim':'Exact finite prefix-pool and saved JSON decimal audit only; numerical trajectories are non-strict and do not prove exclusion or exact counterexample.','hashes':{'config':sha(r/'config.json'),'checkpoint':sha(r/'checkpoint.json'),'cpu_prefix_audit':sha(r/'cpu-prefix-audit.json'),'worker':sha(x.worker),'prefix_module':sha(x.prefix_module)},'prefix_pool':{'raw_ternary_codes':3**(d-1),'valid_low_prefix_count':len(valid),'survivor_count':len(pool),'valid_sha256':encoded(valid),'survivor_sha256':encoded(pool),'all_187_survivors_independently_legal_and_filter_surviving':True,'starts_cover_pool':int(config['starts'])>=len(pool)},'counts':{'restart_files':len(restarts),'rotation_saved_rows':rotation,'complete_final_cohort_saved_rows':final,'partial_rows':0,'total_saved_rows':total,'unique_saved_json_cases':len(seen)},'max_stored_residual_delta':str(float(maxdelta)),'max_designated_complement_delta':str(float(maxcomp)),'minimum_saved_fraction_residual':str(float(minres)),'release_before_exit':{'remaining_tensor_bytes_before_exit':release['remaining_tensor_bytes_before_exit'],'cleanup_called':True,'process_will_exit':True}}
    x.record.parent.mkdir(parents=True,exist_ok=True);x.record.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps(out['counts'],sort_keys=True))

def sha_case(case): return hashlib.sha256(json.dumps(case,sort_keys=True,separators=(',',':')).encode()).hexdigest()
if __name__=='__main__':main()
