#!/usr/bin/env python3
import hashlib,importlib.util,json,math,time
from pathlib import Path
from collections import Counter
T=Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials');HERE=Path(__file__).resolve().parent
I=T/'round038-macbook-gpu52-direct/inputs.json';E=T/'agents/gpu_partial39/native_production42b/takeover44/parent-evidence-task14-23.json';M=E.with_name('handoff-after-task23.json')
V=T/'agents/section_eigen39/native_independent42/reconstruct_verify.py'
SW='97918305c6c2bcc02db02c36a84bf970264a5d387c2491cf3cae4f7233a3711d';BIN='8235c31db552c331e3552044347d3ba784029e7ab7f027bb2440cff23b0ea457'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def existing(info):
 p=Path(info['path']);assert info['exists'] and p.stat().st_size==info['bytes'] and sha(p)==info['sha256'];return p

def main():
 start=time.monotonic();doc=read(I);man=read(M);assert man['source_inputs_sha256']==sha(I)
 spans={i:[] for i in range(len(doc['tasks']))};layerstats=[]
 for layer in man['layers']:
  p=Path(layer['path']);assert sha(p)==layer['sha256'];count=0
  for line in p.read_text().splitlines():
   if not line.startswith('{'):assert '"event"' not in line;continue
   e=json.loads(line);kind=e.get('event',e.get('kind'))
   if kind=='batch':t=e['task'];lo=e['first'];hi=e['last']+1;assert hi-lo==e['count']
   elif kind=='streamed_batch_numeric_ranking':
    t=e['task'];lo,hi=e['range'];assert e['coverage_credit'] is True and hi-lo==e['count']
   else:continue
   assert all(type(x)is int for x in [t,lo,hi]) and 0<=lo<hi<=doc['tasks'][t]['masks'];spans[t].append((lo,hi));count+=hi-lo
  layerstats.append(dict(name=layer['name'],masks=count,sha256=layer['sha256']))
 for t,rows in spans.items():
  end=0
  for lo,hi in sorted(rows):assert lo==end,(t,lo,end);end=hi
  assert end==man['task_next'][str(t)]==doc['tasks'][t]['masks']
 assert sha(V)=='189dd0feb03e07aa48b69c25fa5ecd6172c762ca83b8434499e6fd1edfd64226'
 spec=importlib.util.spec_from_file_location('independent',V);v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)
 evidence=read(E);assert [r['task'] for r in evidence['rows']]==list(range(14,24));summary=[];saved=[]
 for row in evidence['rows']:
  t=row['task'];original=row['original'];mp=existing(original['manifest.json']);mm=read(mp);ip=existing(original[f'task{t}-input.json']);inp=read(ip)
  assert mm['coverage_credit'] is True and mm['task']==t and mm['range']==[0,doc['tasks'][t]['masks']]
  assert mm['input_sha256']==sha(ip) and mm['source_sha256']==SW and inp['source_inputs_sha256']==sha(I) and inp['groups']==doc['tasks'][t]['groups']
  if row['recovered']:
   lp=existing(row['recovered']['ledger']);rp=existing(row['recovered']['receipt']);jp=existing(row['remote_retrieved']['journal']);fp=existing(row['remote_retrieved']['final'])
  else:
   lp=existing(original['ledger.jsonl']);rp=existing(original['receipt.json']);jp=existing(original['remote-stream.jsonl']);fp=existing(original['remote-final.json'])
  rr=read(rp);ff=read(fp);assert rr['status']=='complete' and rr['coverage_credit'] is True and rr['requested_range']==mm['range'] and rr['task']==t
  assert rr['ledger_sha256']==sha(lp) and rr['remote_journal_sha256']==sha(jp)
  if 'manifest_sha256' in rr:assert rr['manifest_sha256']==sha(mp)
  if 'remote_final_sha256' in rr:assert rr['remote_final_sha256']==sha(fp)
  for k,val in [('input_sha256',sha(ip)),('source_sha256',SW),('binary_sha256',BIN)]:assert rr['remote_'+k]==ff[k]==val
  assert ff['command_status']=='completed' and ff['task']==str(t)
  stream=[json.loads(x) for x in jp.read_text().splitlines()];ledger=[json.loads(x) for x in lp.read_text().splitlines()];assert len(stream)==len(ledger)==len(mm['batches'])==rr['verified_batches']
  counts=Counter();maxerr=0;task=v.task_from_source(doc,t)
  for j,e in enumerate(stream):
   lo,n=mm['batches'][j];f=ledger[j];assert e['event']=='batch' and e['command_status']=='completed' and e['task']==str(t) and e['range_start']==lo and e['range_count']==e['finite_count']==n
   for k in ['groups','ncoef','nfreq']:assert e[k]==mm['dims'][k]
   for k,val in [('input_sha256',sha(ip)),('source_sha256',SW),('binary_sha256',BIN)]:assert e[k]==f[k]==val
   assert f['range']==[lo,lo+n] and f['count']==n and f['coverage_credit'] is True
   ix=e['top_indices'];scores=e['top_scores'];assert len(ix)==len(scores)==len(set(ix))==4 and all(type(x)is int and 0<=x<n for x in ix) and all(math.isfinite(s) for s in scores)
   masks=[lo+x for x in ix];assert masks==[x['mask'] for x in f['candidates_cpu_reconstructed']]
   for mask in masks:
    z=v.reconstruct(task,mask,1e-8,1e-8,1e-5);assert z['finite'];counts[z['classification']]+=1;maxerr=max(maxerr,z['max_product_target_error']);saved.append(dict(task=t,mask=mask,classification=z['classification'],min_coefficient=z['min_coefficient'],product_error=z['max_product_target_error']))
  summary.append(dict(task=t,batches=len(stream),saved_candidates=sum(counts.values()),classifications=dict(counts),max_product_error=maxerr,ledger_sha256=sha(lp),receipt_sha256=sha(rp)))
 out=dict(status='passed_parent_full_range_and_task14_23_numerical_replay',total_received_masks=sum(s['masks'] for s in layerstats),layers=layerstats,tasks=summary,scope='All 24 tasks numerical range union only; independent numerical reconstruction covers saved task14..23 candidates. Unknown remains unknown. Not a certificate or general exclusion.',source_sha256=sha(Path(__file__)),manifest_sha256=sha(M),index_sha256=sha(E),checker_sha256=sha(V),seconds=time.monotonic()-start)
 (HERE/'parent-summary.json').write_text(json.dumps(out,indent=2)+'\n');(HERE/'parent-candidate-diagnostics.json').write_text(json.dumps(saved,separators=(',',':'))+'\n');print(json.dumps(out))
if __name__=='__main__':main()
