#!/usr/bin/env python3
"""Independent per-total semantic verifier for the fresh 47/48 streamed run."""
from __future__ import annotations
import argparse, hashlib, importlib.util, json, time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

HERE=Path(__file__).resolve().parent; AUDITED=HERE/'audited_chunk_verifier.py'; AUDITED_SHA='a557bf120c6bed03adbd8b2ce153a88bf76e678ed94c6bc8815d14dc33edbe4b'; CAP=6
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for x in iter(lambda:f.read(1<<20),b''):h.update(x)
 return h.hexdigest()
def audited():
 if sha(AUDITED)!=AUDITED_SHA:raise RuntimeError('audited bytes changed')
 s=importlib.util.spec_from_file_location('audited',AUDITED);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def rows(p):return [json.loads(x) for x in Path(p).read_text().splitlines() if x]
def key(x):return int(x['total']),int(x['m']),int(x['n']),int(x['index']),tuple(x['I'])
def payload(frozen,remote,jobs,plan_sha):
 text=Path(frozen).read_text();pre,mark,_=text.partition('\ndef main():')
 if not mark:raise RuntimeError('frozen main marker missing')
 return ('SOURCE_SHA='+repr(sha(frozen))+'\n'+pre+'\n'+Path(remote).read_text()+'\nrun_support48(jobs='+json.dumps(jobs,separators=(',',':'))+', plan_sha256='+repr(plan_sha)+', source_sha256=SOURCE_SHA, workers=15, seconds=14400)\n').encode()
def expected(a,total):return {(total,m,total-m,i,s) for m in range(1,total//2+1) for i,s in enumerate(a.independent_supports(m,CAP))}
def manifest(run,plan_path,frozen,remote):
 a=audited(); plan=json.loads(Path(plan_path).read_text()); ev=rows(Path(run)/'events.jsonl'); errors=[]; plan_sha=sha(plan_path)
 if plan.get('format')!='a30-total47-48-cap6-v1' or plan.get('parameters')!={'cap':6,'remote_seconds_limit':14400,'remote_workers':15,'totals':[47,48]}:errors.append('plan_contract')
 if plan.get('frozen_support_sha256')!=sha(frozen) or plan.get('job_total')!=100005:errors.append('plan_frozen_or_size')
 want=expected(a,47)|expected(a,48); jobs={key(x):x for x in plan.get('jobs',[])}
 if len(jobs)!=len(plan.get('jobs',[])) or set(jobs)!=want:errors.append('plan_not_full_independent_47_48_table')
 launch=json.loads((Path(run)/'launch-manifest.json').read_text()); body=payload(frozen,remote,plan.get('jobs',[]),plan_sha)
 if any(launch.get(k)!=v for k,v in {'format':'a30-support48-launch-v1','plan_sha256':plan_sha,'frozen_support_sha256':sha(frozen),'remote_sha256':sha(remote),'payload_sha256':hashlib.sha256(body).hexdigest(),'selected_jobs':100005,'workers':15,'seconds_limit':14400}.items()):errors.append('launch_binding')
 if (Path(run)/'payload.py').read_bytes()!=body:errors.append('payload_bytes')
 receipt=json.loads((Path(run)/'exit-receipt.json').read_text())
 if not(receipt.get('returncode')==0 and receipt.get('selected_set_complete') is True and receipt.get('chunks_persisted')==100005 and receipt.get('terminal_event')=='support48_selected_complete'):errors.append('receipt')
 if len(ev)<4 or ev[0].get('event')!='support48_local_started' or ev[0].get('plan_sha256')!=plan_sha:errors.append('local_start')
 starts=[x for x in ev if x.get('event')=='support48_remote_started']; workers=[x for x in ev if x.get('event')=='support48_remote_workers']
 if len(starts)!=1 or starts[0].get('frozen_source_sha256')!=sha(frozen) or starts[0].get('workers')!=15 or len(workers)!=1 or len(workers[0].get('worker_pids',[]))!=15:errors.append('remote_start_workers')
 chunks=[x for x in ev if x.get('event')=='chunk']; seen=set()
 for x in chunks:
  try:k=key(x)
  except Exception:errors.append('malformed_chunk');continue
  if k not in want:errors.append('outside_table');continue
  if k in seen:errors.append('duplicate_chunk');continue
  seen.add(k)
 if seen!=want:errors.append('events_not_full_exact_coverage')
 term=ev[-1] if ev else {}
 if term.get('event')!='support48_selected_complete' or term.get('workers_joined') is not True or term.get('submitted')!=100005 or term.get('finished')!=100005:errors.append('terminal')
 return {'errors':errors,'plan_sha256':plan_sha,'events_sha256':sha(Path(run)/'events.jsonl'),'frozen_sha256':sha(frozen),'remote_sha256':sha(remote),'payload_sha256':hashlib.sha256(body).hexdigest(),'receipt_sha256':sha(Path(run)/'exit-receipt.json'),'audited_sha256':sha(AUDITED),'chunks':chunks,'plan':plan}
def worker(run,event):
 a=audited();e=[];r=a.verify_chunk(Path(run)/f"t{event['total']:02d}",event,CAP,e)
 return {'job':r['job'],'errors':e,'verified_rows':r['verified_rows'],'counts':dict(r['counts']),'gzip_sha256':event['gzip_sha256'],'raw_sha256':event['raw_sha256']}
def ckey(total,event):return f"t{total:02d}-{event['m']:02d}-{event['n']:02d}-{event['index']:06d}.json"
def cache_ok(path,run,event,binding):
 try:
  x=json.loads(path.read_text())
  if x['binding']!=binding or x['event']!={'gzip_sha256':event['gzip_sha256'],'raw_sha256':event['raw_sha256']}:return None
  target=Path(run)/f"t{event['total']:02d}"/f"m{event['m']:02d}n{event['n']:02d}"/f"i{event['index']:06d}.jsonl.gz"
  return x['result'] if sha(target)==event['gzip_sha256'] and not x['result']['errors'] else None
 except Exception:return None
def atom(p,x):
 q=p.with_suffix('.tmp');q.write_text(json.dumps(x,sort_keys=True)+'\n');q.replace(p)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('total',type=int,choices=(47,48));ap.add_argument('--run',type=Path,required=True);ap.add_argument('--plan',type=Path,required=True);ap.add_argument('--frozen',type=Path,required=True);ap.add_argument('--remote',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--cache',type=Path,required=True);ap.add_argument('--workers',type=int,default=4);z=ap.parse_args()
 if not Path('/Volumes/data').is_mount() or not 1<=z.workers<=4:raise SystemExit('mount/worker failure')
 began=time.monotonic();m=manifest(z.run,z.plan,z.frozen,z.remote);report={'format':'a30-independent-support48-per-total-v1','total':z.total,'workers':z.workers,'manifest':{k:v for k,v in m.items() if k not in ('chunks','plan')},'status':'rejected_manifest','chunk_results':[],'counts':{}}
 if not m['errors']:
  binding={k:m[k] for k in ('plan_sha256','events_sha256','frozen_sha256','remote_sha256','payload_sha256','receipt_sha256','audited_sha256')}; z.cache.mkdir(parents=True,exist_ok=True); todo=[];tot=Counter(raw=0,canonical=0,excluded=0,fair=0,unknown=0);reuse=0
  selected=[x for x in m['chunks'] if int(x['total'])==z.total]
  for x in selected:
   c=cache_ok(z.cache/ckey(z.total,x),z.run,x,binding)
   if c is None:todo.append(x)
   else:report['chunk_results'].append(c);tot.update(c['counts']);reuse+=1
  with ProcessPoolExecutor(max_workers=z.workers) as pool:
   fs={pool.submit(worker,str(z.run),x):x for x in todo}
   for f in as_completed(fs):
    x=fs[f];r=f.result();report['chunk_results'].append(r);tot.update(r['counts']);atom(z.cache/ckey(z.total,x),{'binding':binding,'event':{'gzip_sha256':x['gzip_sha256'],'raw_sha256':x['raw_sha256']},'result':r})
  report['chunk_results'].sort(key=lambda x:x['job']);report['counts']=dict(tot);report['cache']={'binding':binding,'reused':reuse,'scheduled':len(todo)};report['status']='passed_complete_total_'+str(z.total) if not any(x['errors'] for x in report['chunk_results']) else 'rejected_chunk_content'
 report['elapsed_seconds']=time.monotonic()-began;report['scope']='Fixed total '+str(z.total)+' only; unknown remains unknown; no general conclusion.';z.out.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');print(json.dumps({'status':report['status'],'chunks':len(report['chunk_results']),'counts':report['counts'],'seconds':report['elapsed_seconds']},sort_keys=True));raise SystemExit(0 if report['status']=='passed_complete_total_'+str(z.total) else 2)
if __name__=='__main__':main()
