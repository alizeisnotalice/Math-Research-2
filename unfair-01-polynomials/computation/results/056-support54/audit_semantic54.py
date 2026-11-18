"""Independent local audit of full returned total54 receipt coverage and provenance."""
import hashlib,itertools,json,os,time
from pathlib import Path
from collections import Counter
B=Path(__file__).resolve().parent;T=B.parent;R=B/'support54-remote'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 if not __debug__:raise RuntimeError('no -O')
 out=B/'semantic54-audit.json';assert not out.exists();start=time.monotonic();inp=R/'input';bind={p.name:sha(p) for p in inp.iterdir() if p.is_file()};bind.update({'runner':sha(T/'round053-parent/support54_runner.py'),'source':sha(inp/'frozen_support.py'),'verifier':sha(R/'code/audited_chunk_verifier.py'),'total':54,'cap':6})
 assert bind['source']=='d68f1010c61ad1d911749742fef8c790c0a8064bf34512b402e7e65e4703e710' and bind['verifier']=='a557bf120c6bed03adbd8b2ce153a88bf76e678ed94c6bc8815d14dc33edbe4b'
 progress=json.loads((R/'progress54.json').read_text());end=json.loads((R/'exit54.json').read_text());launch=json.loads((R/'launch54.json').read_text())
 assert progress['status']=='complete' and progress['expected']==progress['verified']==101583 and progress['remaining']==0 and progress['failures']==[] and progress['binding']==bind
 assert end['returncode']==0 and end['child_waited'] and end['remaining']==0 and end['status']=='complete' and end['verified']==101583 and end['pid']==launch['pid']==progress['pid'] and launch['runner_sha256']==bind['runner']
 expected={}
 for m in range(1,28):
  i=0
  for size in range(2,min(m+1,6)+1):
   for inner in itertools.combinations(range(1,m),size-2):expected[m,i]=(0,*inner,m);i+=1
 plan=json.loads((inp/'plan.json').read_text());assert plan['job_total']==len(plan['jobs'])==185264
 actual={}
 for j in plan['jobs']:
  if j['total']!=54:continue
  k=j['m'],j['index'];assert k not in actual and j['n']==54-j['m'];actual[k]=tuple(j['I'])
 assert actual==expected and len(expected)==101583
 old=json.loads((inp/'old-exit.json').read_text());new=json.loads((inp/'new-exit.json').read_text());ol=json.loads((inp/'old-launch.json').read_text())
 assert old['returncode']==255 and not old['selected_set_complete'] and old['chunks_persisted']==94465
 assert new['status']=='complete' and not new['failures'] and new['remaining']==0 and new['persisted_chunks_total']==90799 and new['submitted']==new['finished']==90798
 assert ol['plan_sha256']==sha(inp/'plan.json') and ol['payload_sha256']==sha(inp/'old-payload.py') and ol['frozen_support_sha256']==bind['source']
 events={};origin=Counter();newkeys=set()
 for label in ('old','new'):
  for line in (inp/f'{label}-events.jsonl').open():
   e=json.loads(line)
   if e.get('event')!='chunk' or e.get('total')!=54:continue
   k=e['m'],e['index'];assert k not in events and tuple(e['I'])==expected[k] and e['n']==54-k[0];events[k]=e;origin[label]+=1
   if label=='new':newkeys.add(k)
 assert origin=={'old':10784,'new':90799} and events.keys()==expected.keys()
 comp=json.loads((inp/'complement.json').read_text());assert comp['missing_total54_jobs']==90799 and {(j['m'],j['index']):tuple(j['I']) for j in comp['missing_jobs']}=={k:expected[k] for k in newkeys}
 counts=Counter();seen=set()
 for p in (R/'cache/t54').glob('m*/i*.json'):
  k=int(p.parent.name[1:]),int(p.stem[1:]);assert k in expected and k not in seen;seen.add(k);x=json.loads(p.read_text());e=events[k]
  assert x['status']=='passed' and x['errors']==[] and x['binding']==bind and x['event_sha256']==hashlib.sha256(json.dumps(e,sort_keys=True).encode()).hexdigest()
  z=x['result'];assert z['counts']==e['counts'] and z['verified_rows']==e['counts']['canonical'] and z['raw_support_count']==e['counts']['raw'] and z['job']==f'm{k[0]:02}n{54-k[0]:02}/i{k[1]:06}'
  counts.update(e['counts'])
 assert seen==expected.keys() and counts['unknown']==0
 raw=B/'raw54-audit/summary.json';assert raw.exists(),'raw coverage audit must finish first';a=json.loads(raw.read_text());assert a['status']=='passed' and a['chunks']==101583 and a['counts']==dict(counts)
 union={}
 for line in (B/'raw54-input/union-events.jsonl').open():
  e=json.loads(line);e.pop('_run');k=e['m'],e['index'];assert k not in union;union[k]=e
 assert union==events
 result={'status':'passed','chunks':len(seen),'counts':dict(counts),'binding':bind,'origin':dict(origin),'raw_audit_sha256':sha(raw),'remote_exit_sha256':sha(R/'exit54.json'),'auditor_sha256':sha(Path(__file__)),'seconds':time.monotonic()-start,'scope':'Complete total54 semantic receipt coverage plus independent local full gzip/raw byte coverage; no claim beyond cap6 total54.'}
 out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':main()
