#!/usr/bin/env python3
"""Audit saved decimal samples; the previously certified prefix pool is reused."""
import argparse,hashlib,importlib.util,json
from fractions import Fraction
from pathlib import Path

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);ap.add_argument('--verifier',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();r=a.run
 s=importlib.util.spec_from_file_location('independent',a.verifier);v=importlib.util.module_from_spec(s);s.loader.exec_module(v)
 sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
 c=json.loads((r/'config.json').read_text());p=json.loads((r/'prefix-plan.json').read_text());x=json.loads((r/'result.json').read_text())
 assert sha(r/'worker.py')==c['script_sha256'] and sha(r/'prefix-plan.json')==c['prefix_plan_sha256']
 assert x['complete'] and x['last_completed_step']==6000
 pool=set(p['survivor_codes']);assert len(pool)==46011
 seen={p['survivor_codes'][(epoch*8192+i)%len(pool)] for epoch in range(6) for i in range(8192)};assert seen==pool
 ans=[]
 for case in x['cases']:
  res,delta,comp=v.validate_case(case,64,104,15,Fraction(1,1000),pool);assert res>0
  ans.append({'exact_decimal_residual':str(res),'float_residual':float(res),'complement_error':str(comp)})
 y={'status':'passed_16_saved_decimal_samples_and_window_assignment','scope':'Exact decimal sample audit; no exact counterexample and no exclusion. Pool completeness inherited from previously verified d15 input, not rerun here.','sample_count':len(ans),'coverage_prefixes':len(seen),'residuals':ans,'bindings':{n:sha(r/n) for n in ['worker.py','config.json','prefix-plan.json','result.json','release.json']},'verifier_sha256':sha(a.verifier)}
 a.out.write_text(json.dumps(y,indent=2)+'\n');print(json.dumps({'samples':len(ans),'prefixes':len(seen),'min_residual':min(v['float_residual'] for v in ans)}))
if __name__=='__main__':main()
