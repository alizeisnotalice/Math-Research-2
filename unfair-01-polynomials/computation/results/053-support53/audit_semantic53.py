"""Audit complete fetched total53 semantic receipts and original interrupted provenance."""
import argparse,hashlib,itertools,json,os,time
from collections import Counter
from pathlib import Path

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--temp',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();T=a.temp;B=Path(__file__).parent;local=B/'support53-remote';run=T/'round043-macbook-support54-corrected/run';planpath=T/'agents/support54_macbook43/plan/support54-plan.json';source=planpath.parent/'frozen_support.py';runner=T/'round050-parent/support53/support53_batch_runner.py';verifier=T/'round050-parent/support53/code__audited_chunk_verifier.py';start=time.monotonic()
 assert not a.out.exists(),'fresh output required'
 plan=json.loads(planpath.read_text());assert len(plan['jobs'])==185264
 expected={}
 for m in range(1,27):
  i=0
  for size in range(2,min(m+1,6)+1):
   for inner in itertools.combinations(range(1,m),size-2):expected[m,i]=(0,*inner,m);i+=1
 planned={}
 for j in plan['jobs']:
  if j['total']!=53:continue
  key=(j['m'],j['index']);assert key not in planned and j['n']==53-j['m'];planned[key]=tuple(j['I'])
 assert expected==planned and len(expected)==83681
 events={}
 for l in (run/'events.jsonl').open():
  e=json.loads(l)
  if e.get('event')!='chunk' or e.get('total')!=53:continue
  key=e['m'],e['index'];assert key not in events and tuple(e['I'])==expected[key] and e['n']==53-e['m'];events[key]=e
 assert events.keys()==expected.keys()
 original=json.loads((run/'exit-receipt.json').read_text());assert original['returncode']==255 and not original['selected_set_complete'] and original['chunks_persisted']==94465
 launch=json.loads((run/'launch-manifest.json').read_text());assert launch['plan_sha256']==sha(planpath) and launch['payload_sha256']==sha(run/'payload.py') and launch['frozen_support_sha256']==sha(source)
 binding={'cap':6,'events':sha(run/'events.jsonl'),'launch':sha(run/'launch-manifest.json'),'payload':sha(run/'payload.py'),'plan':sha(planpath),'runner':sha(runner),'source':sha(source),'terminal':sha(run/'exit-receipt.json'),'total':53,'verifier':sha(verifier)}
 progress=json.loads((local/'progress53.json').read_text());end=json.loads((local/'exit53.json').read_text());launch53=json.loads((local/'launch53.json').read_text())
 assert progress['status']=='complete' and progress['remaining']==0 and progress['verified']==progress['expected']==83681 and not progress['failures'] and progress['binding']==binding
 assert end['returncode']==0 and end['child_waited'] is True and end['status']=='complete' and end['verified']==83681 and end['remaining']==0 and end['pid']==launch53['pid']==progress['pid']
 assert launch53['runner_sha256']==sha(runner)
 ct=Counter();count=0;cache=local/'cache/t53';actual=set()
 for p in cache.glob('m*/i*.json'):
  key=(int(p.parent.name[1:]),int(p.stem[1:]));assert key not in actual and key in expected;actual.add(key);x=json.loads(p.read_text());e=events[key]
  assert x['status']=='passed' and x['errors']==[] and x['binding']==binding
  assert x['event_sha256']==hashlib.sha256(json.dumps(e,sort_keys=True).encode()).hexdigest()
  assert x['result']['counts']==e['counts'] and x['result']['verified_rows']==e['counts']['canonical'] and x['result']['raw_support_count']==e['counts']['raw']
  assert x['result']['job']==f'm{key[0]:02}n{53-key[0]:02}/i{key[1]:06}'
  ct.update(e['counts']);count+=1
 assert actual==expected.keys() and ct['unknown']==0
 raw=json.loads((B/'raw53-audit/summary.json').read_text());assert raw['status']=='passed' and raw['counts']==dict(ct) and raw['chunks']==83681 and raw['binding']['events']==binding['events']
 o={'status':'passed','chunks':count,'counts':dict(ct),'binding':binding,'raw_audit_sha256':sha(B/'raw53-audit/summary.json'),'remote_exit_sha256':sha(local/'exit53.json'),'auditor_sha256':sha(Path(__file__)),'seconds':time.monotonic()-start,'scope':'Complete total53 independent semantic cache coverage plus all local raw-byte checks. Original joint53/54 run remains interrupted; no54 semantic claim.'}
 t=a.out.with_suffix('.partial');t.write_text(json.dumps(o,indent=2)+'\n');os.replace(t,a.out);print(json.dumps(o))
if __name__=='__main__':main()
