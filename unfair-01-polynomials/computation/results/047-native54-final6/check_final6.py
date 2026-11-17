#!/usr/bin/env python3
"""Bind received native streams and independently reconstruct saved candidates."""
import hashlib,importlib.util,json,math,platform,time
from pathlib import Path
from collections import Counter
T=Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials');HERE=Path(__file__).resolve().parent
D=T/'agents/gpu_partial39/degree54_native45'
I=D/'inputs-degree54.json';E=D/'complete-index-task0-23.json'
V=T/'agents/section_eigen39/native_independent42/reconstruct_verify.py'
SW='fa197902d6ffcdc19baa1bb80dad568c5f06750fe0a4a84f978cb453b864e913'
BIN='e9c2f3677eca4383381cda25c9ffceea8324578a1096eb1decb04e9c313f6445'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())

def main():
 start=time.monotonic();doc=read(I);index=read(E)
 assert doc['degree']==54 and len(doc['tasks'])==24
 assert sha(index['source_path'])==index['source_sha256']==SW
 assert sha(V)=='189dd0feb03e07aa48b69c25fa5ecd6172c762ca83b8434499e6fd1edfd64226'
 spec=importlib.util.spec_from_file_location('independent',V);v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)
 assert [r['task'] for r in index['received']]==list(range(24))
 summary=[];saved=[];total_masks=0
 for row in index['received']:
  if row['task']<18:continue
  t=row['task'];paths={k:Path(row[k+'_path']) for k in row['sha256']}
  for k,p in paths.items():assert sha(p)==row['sha256'][k],(t,k)
  mm=read(paths['manifest']);inp=read(paths['input']);rr=read(paths['receipt']);ff=read(paths['remote_final']);att=read(paths['launch_attestation'])
  assert mm['coverage_credit'] is True and mm['task']==t and mm['range']==[0,doc['tasks'][t]['masks']]
  assert mm['input_sha256']==sha(paths['input']) and mm['source_sha256']==SW
  assert inp['source_inputs_sha256']==sha(I) and inp['groups']==doc['tasks'][t]['groups'] and inp['task']==t
  assert att['coverage_credit'] is True and att['range']==mm['range'] and att['task']==t
  assert att['manifest_sha256']==sha(paths['manifest']) and att['input_sha256']==sha(paths['input'])
  assert att['source_sha256']==SW and att['expected_binary_sha256']==BIN
  assert rr['status']=='complete' and rr['coverage_credit'] is True and rr['requested_range']==mm['range'] and rr['task']==t
  assert rr['ledger_sha256']==sha(paths['ledger']) and rr['remote_journal_sha256']==sha(paths['remote_journal'])
  for k in ['manifest','remote_final','launch_attestation']:
   if k+'_sha256' in rr:assert rr[k+'_sha256']==sha(paths[k])
  for k,val in [('input_sha256',sha(paths['input'])),('source_sha256',SW),('binary_sha256',BIN)]:assert rr['remote_'+k]==ff[k]==val
  assert ff['command_status']=='completed' and ff['task']==str(t)
  stream=[json.loads(x) for x in paths['remote_journal'].read_text().splitlines()]
  ledger=[json.loads(x) for x in paths['ledger'].read_text().splitlines()]
  assert len(stream)==len(ledger)==len(mm['batches'])==rr['verified_batches']==row['event_count']
  assert len(ff['batches'])==len(stream)
  assert all(dict(b,event='batch')==s for b,s in zip(ff['batches'],stream))
  counts=Counter();maxerr=0;task=v.task_from_source(doc,t);end=0
  for j,e in enumerate(stream):
   lo,n=mm['batches'][j];f=ledger[j];assert lo==end and n>0;end=lo+n
   assert e['event']=='batch' and e['command_status']=='completed' and e['task']==str(t)
   assert e['range_start']==lo and e['range_count']==e['finite_count']==n
   for k in ['groups','ncoef','nfreq']:assert e[k]==mm['dims'][k]
   for k,val in [('input_sha256',sha(paths['input'])),('source_sha256',SW),('binary_sha256',BIN)]:assert e[k]==f[k]==val
   assert f['range']==row['ranges'][j]==[lo,lo+n] and f['count']==n and f['coverage_credit'] is True
   ix=e['top_indices'];scores=e['top_scores']
   assert len(ix)==len(scores)==len(set(ix))==4 and all(type(x)is int and 0<=x<n for x in ix) and all(math.isfinite(s) for s in scores)
   masks=[lo+x for x in ix];assert masks==[x['mask'] for x in f['candidates_cpu_reconstructed']]
   for mask in masks:
    z=v.reconstruct(task,mask,1e-8,1e-8,1e-5);assert z['finite'];counts[z['classification']]+=1;maxerr=max(maxerr,z['max_product_target_error'])
    saved.append(dict(task=t,mask=mask,classification=z['classification'],min_coefficient=z['min_coefficient'],product_error=z['max_product_target_error']))
  assert end==doc['tasks'][t]['masks'];total_masks+=end
  summary.append(dict(task=t,masks=end,batches=len(stream),saved_candidates=sum(counts.values()),classifications=dict(counts),max_product_error=maxerr))
 out=dict(status='passed_parent_task18_23_binding_and_numerical_replay',python=platform.python_version(),masks=total_masks,tasks=summary,scope='Only degree54 tasks18..23 replayed here; 0..17 were certified numerically in round45. Final objects are matched to journal and input/source/binary; task19..23 receipts do not directly contain final-file hashes. Numerical reconstruction is not an exact certificate or exclusion.',source_sha256=sha(Path(__file__)),index_sha256=sha(E),inputs_sha256=sha(I),checker_sha256=sha(V),seconds=time.monotonic()-start)
 (HERE/'parent-summary.json').write_text(json.dumps(out,indent=2)+'\n')
 (HERE/'candidate-diagnostics.json').write_text(json.dumps(saved,separators=(',',':'))+'\n')
 print(json.dumps(out))
if __name__=='__main__':main()
