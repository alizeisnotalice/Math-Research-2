#!/usr/bin/env python3
"""Independent 49/50 old+resume partition manifest and semantic worker adapter."""
import argparse,hashlib,importlib.util,json
from pathlib import Path
HERE=Path(__file__).resolve().parent; AUD=HERE/'audited_chunk_verifier.py'; PIN='a557bf120c6bed03adbd8b2ce153a88bf76e678ed94c6bc8815d14dc33edbe4b'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load():
 if sha(AUD)!=PIN:raise RuntimeError('audited verifier changed')
 s=importlib.util.spec_from_file_location('aud',AUD);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def k(x):return int(x['total']),int(x['m']),int(x['n']),int(x['index']),tuple(x['I'])
def events(p):return [json.loads(x) for x in Path(p).read_text().splitlines() if x]
def check(old,resume,plan,resume_plan):
 a=load(); P=json.loads(Path(plan).read_text());R=json.loads(Path(resume_plan).read_text());err=[]
 want={(t,m,t-m,i,s) for t in (49,50) for m in range(1,t//2+1) for i,s in enumerate(a.independent_supports(m,6))}
 if {k(x) for x in P['jobs']}!=want or len(P['jobs'])!=123859:err.append('fresh_plan_not_independent_full_table')
 if len({k(x) for x in R['jobs']})!=len(R['jobs']):err.append('duplicate_resume_plan_jobs')
 oe=[x for x in events(Path(old)/'events.jsonl') if x.get('event')=='chunk']; ne=[x for x in events(Path(resume)/'events.jsonl') if x.get('event')=='chunk']
 os,ns={k(x) for x in oe},{k(x) for x in ne}; missing={k(x) for x in R['jobs']}
 if len(os)!=len(oe) or len(ns)!=len(ne):err.append('duplicate_event_chunk')
 if len(os)!=58036 or len(ns)!=65823:err.append('unexpected_side_counts')
 if ns!=missing or os&ns or os|ns!=want:err.append('old_resume_not_exact_complement')
 receipt=json.loads((Path(resume)/'exit.json').read_text())
 if not(receipt.get('returncode')==0 and receipt.get('selected_set_complete') and receipt.get('chunks_persisted')==65823 and receipt.get('all_owned_processes_reaped')):err.append('resume_exit_contract')
 return {'errors':err,'old_chunks':len(oe),'resume_chunks':len(ne),'expected':len(want),'old_events_sha256':sha(Path(old)/'events.jsonl'),'resume_events_sha256':sha(Path(resume)/'events.jsonl'),'old_payload_sha256':sha(Path(old)/'payload.py'),'old_receipt_sha256':sha(Path(old)/'exit-receipt.json'),'resume_payload_sha256':sha(Path(resume)/'payload.py'),'resume_receipt_sha256':sha(Path(resume)/'exit.json'),'fresh_frozen_sha256':sha(Path(plan).parent/'frozen_support.py'),'fresh_remote_sha256':sha('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/support50_macbook41/remote_support50.py'),'resume_frozen_sha256':sha('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/gpu_partial39/support50_resume41/frozen_support.py'),'resume_remote_sha256':sha('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/gpu_partial39/support50_resume41/remote_support50.py'),'fresh_plan_sha256':sha(plan),'resume_plan_sha256':sha(resume_plan),'audited_sha256':sha(AUD)}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--old',type=Path,required=True);ap.add_argument('--resume',type=Path,required=True);ap.add_argument('--plan',type=Path,required=True);ap.add_argument('--resume-plan',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);z=ap.parse_args();r=check(z.old,z.resume,z.plan,z.resume_plan);z.out.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n');print(json.dumps(r,sort_keys=True));raise SystemExit(2 if r['errors'] else 0)
if __name__=='__main__':main()
