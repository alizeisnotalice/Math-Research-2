#!/usr/bin/env python3
"""Independent total53 semantic replay from interrupted53/54 production; exact complete53 slice only."""
import argparse,concurrent.futures as cf,hashlib,importlib.util,itertools,json,multiprocessing,os,sys,time
from pathlib import Path
if not __debug__:raise RuntimeError('Optimized Python disables required checks')
PIN='a557bf120c6bed03adbd8b2ce153a88bf76e678ed94c6bc8815d14dc33edbe4b'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def eh(e):return hashlib.sha256(json.dumps(e,sort_keys=True).encode()).hexdigest()
def atom(p,o):
 p.parent.mkdir(parents=True,exist_ok=True);q=p.with_suffix('.partial')
 with q.open('w') as f:f.write(json.dumps(o,sort_keys=True)+'\n');f.flush();os.fsync(f.fileno())
 q.replace(p)
def worker(root,e,verifier,binding,cache):
 s=importlib.util.spec_from_file_location('aud',verifier);v=importlib.util.module_from_spec(s);s.loader.exec_module(v)
 errors=[];r=v.verify_chunk(root/f't{e["total"]}',e,6,errors)
 o={'status':'passed' if not errors else 'rejected','binding':binding,'event_sha256':eh(e),'errors':errors,'result':r}
 atom(cache,o);return o

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--total',type=int,choices=(53,),required=True);ap.add_argument('--workers',type=int,default=13);ap.add_argument('--seconds',type=int,default=14400);ap.add_argument('--limit',type=int);a=ap.parse_args();start=time.monotonic()
 root=a.root;planpath=root/'plan/support54-plan.json';vp=root/'code/audited_chunk_verifier.py';assert sha(vp)==PIN
 plan=json.loads(planpath.read_text());assert plan['format']=='a30-total53-54-cap6-v1' and plan['parameters']=={'cap':6,'remote_seconds_limit':14400,'remote_workers':15,'totals':[53,54]} and plan['job_total']==len(plan['jobs'])==185264
 assert len({(j['total'],j['m'],j['index']) for j in plan['jobs']})==185264
 assert all(j['total'] in (53,54) for j in plan['jobs'])
 launch=json.loads((root/'run/launch-manifest.json').read_text());terminal=json.loads((root/'run/exit-receipt.json').read_text())
 assert launch['plan_sha256']==sha(planpath) and launch['payload_sha256']==sha(root/'run/payload.py') and launch['frozen_support_sha256']==plan['frozen_support_sha256']
 # Original 53/54 production was interrupted; only the complete 53 slice is replayed.
 assert terminal['returncode']==255 and terminal['selected_set_complete'] is False and terminal['chunks_persisted']==94465
 assert a.total==53
 assert sha(root/'plan/frozen_support.py')==plan['frozen_support_sha256']
 expected={}
 for m in range(1,a.total//2+1):
  ix=0
  for size in range(2,min(m+1,6)+1):
   for inside in itertools.combinations(range(1,m),size-2):expected[(m,ix)]=(0,*inside,m);ix+=1
 actual={}
 for j in plan['jobs']:
  if j['total']!=a.total:continue
  k=(j['m'],j['index']);assert k not in actual and j['n']==a.total-j['m'];actual[k]=tuple(j['I'])
 assert actual==expected
 del plan,actual
 events={};terminal_events=[]
 with (root/'run/events.jsonl').open() as f:
  for line in f:
   e=json.loads(line)
   if e.get('event')=='support54_selected_complete':terminal_events.append(e)
   if e.get('event')!='chunk' or e.get('total')!=a.total:continue
   k=(e['m'],e['index']);assert k in expected and k not in events and e['n']==a.total-e['m'] and tuple(e['I'])==expected[k];events[k]=e
 assert events.keys()==expected.keys()
 assert not terminal_events and len(expected)==83681
 # No overall production completion is claimed. Exact layer-53 support/event coverage above is required.
 bind={'plan':sha(planpath),'events':sha(root/'run/events.jsonl'),'source':sha(root/'plan/frozen_support.py'),'verifier':PIN,'runner':sha(Path(__file__)),'total':a.total,'cap':6,'payload':sha(root/'run/payload.py'),'launch':sha(root/'run/launch-manifest.json'),'terminal':sha(root/'run/exit-receipt.json')}
 def paths(k):
  m,i=k;return root/f'run/t{a.total}/m{m:02}n{a.total-m:02}/i{i:06}.jsonl.gz',root/f'cache/t{a.total}/m{m:02}/i{i:06}.json'
 pending=set(events);done=0;failed=[];cached=0
 for k in sorted(pending):
  p,c=paths(k)
  if not c.is_file():continue
  try:
   r=json.loads(c.read_text());e=events[k]
   assert r['status']=='passed' and not r['errors'] and r['binding']==bind and r['event_sha256']==eh(e) and r['result']['counts']==e['counts'] and sha(p)==e['gzip_sha256']
   pending.remove(k);done+=1;cached+=1
  except (OSError,KeyError,ValueError,AssertionError):pass
 def progress(status):
  o={'status':status,'total':a.total,'expected':len(events),'verified':done,'cached':cached,'remaining':len(events)-done,'failures':failed,'binding':bind,'pid':os.getpid(),'python':sys.version,'workers':a.workers,'elapsed':time.monotonic()-start}
  atom(root/f'progress{a.total}.json',o);print(json.dumps({k:o[k] for k in ('status','total','verified','remaining','elapsed')}),flush=True)
 progress('running');submitted=0;last=0
 with cf.ProcessPoolExecutor(max_workers=a.workers,mp_context=multiprocessing.get_context('fork')) as pool:
  active={}
  while pending or active:
   expired=time.monotonic()-start>=a.seconds or (a.limit is not None and submitted>=a.limit)
   if not expired and not failed:
    for k in sorted(pending):
     if len(active)>=2*a.workers or (a.limit is not None and submitted>=a.limit):break
     p,c=paths(k)
     if not p.is_file():continue
     # rsync uses a temporary filename and atomic rename. Verifier rechecks both hashes.
     fut=pool.submit(worker,root/'run',events[k],vp,bind,c);active[fut]=k;pending.remove(k);submitted+=1
   if not active:
    if expired or failed:break
    time.sleep(2)
   else:
    ready,_=cf.wait(active,timeout=2,return_when=cf.FIRST_COMPLETED)
    for fut in ready:
     k=active.pop(fut)
     try:
      r=fut.result()
      if r['status']=='passed':done+=1
      else:failed.append({'key':k,'errors':r['errors']})
     except Exception as e:failed.append({'key':k,'error':repr(e)})
   if time.monotonic()-last>=10:progress('running' if not failed else 'draining_after_failure');last=time.monotonic()
 progress('complete' if done==len(events) and not failed else 'rejected' if failed else 'partial')
 if failed:raise SystemExit(2)
if __name__=='__main__':main()
