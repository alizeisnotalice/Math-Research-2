#!/usr/bin/env python3
"""Fail-closed controller for fresh l=4,5 row-orbit production or tiny smoke.

SIGINT/SIGTERM only stops new dispatch.  Already spawned children are waited
and reaped; this program never sends a signal to a child.
"""
from __future__ import annotations
import argparse, concurrent.futures as cf, hashlib, json, os, signal, subprocess, sys, tempfile, time, uuid
from pathlib import Path

N=11; STOP=False; STOP_REASON=None
class Bad(RuntimeError): pass
def now(): return time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def put(p,x):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True); q=p.with_name('.'+p.name+'.tmp')
 with q.open('w') as f:json.dump(x,f,sort_keys=True);f.write('\n')
 os.replace(q,p)
def get(p): return json.loads(Path(p).read_text())
def comb(n,k):
 import math
 return math.comb(n,k)
def total(l): return comb((1<<(N-l))+l-1,l)
def on_signal(sig,_frame):
 global STOP,STOP_REASON
 STOP=True; STOP_REASON=signal.Signals(sig).name
def valid_result(p,l,a,b,source,sid):
 d=get(p)
 need={'status':'completed','mode':'scan-roworbit-j7','N':N,'left_size':l,'orbit_start':a,'orbit_stop_exclusive':b,'next_orbit_rank':b,'orbits_scanned':b-a,'max_abs_Rj7':0,'source_sha256':source,'shard_id':sid}
 if any(d.get(k)!=v for k,v in need.items()): raise Bad(f'checkpoint semantic failure {p}')
 if not isinstance(d.get('weighted_labeled_graphs'),int) or d['weighted_labeled_graphs']<0: raise Bad(f'bad weight {p}')
def task(binary,source,out,session,l,a,b):
 sid=f'rowprod-{session}-l{l}-{a:012d}-{b:012d}'; logs=out/'logs'; results=out/'results'; cmd=['/usr/bin/nice','-n','10',str(binary),'--mode','scan-roworbit-j7','--left-size',str(l),'--start',str(a),'--stop',str(b),'--checkpoint-every','256','--shard-id',sid,'--out-dir',str(results),'--source-sha',source]
 spawn=logs/(sid+'.spawn.json'); wait=logs/(sid+'.wait.json'); reap=logs/(sid+'.reap.json'); error=logs/(sid+'.receipt-error.json')
 begun=now(); child=None; waited=False; rc=None
 try:
  with (logs/(sid+'.stdout')).open('w') as so,(logs/(sid+'.stderr')).open('w') as se:
   child=subprocess.Popen(cmd,stdout=so,stderr=se); put(spawn,{'event':'spawn','utc':begun,'pid':child.pid,'command':cmd})
   rc=child.wait(); waited=True; put(wait,{'event':'wait','utc':now(),'pid':child.pid,'returncode':rc}); put(reap,{'event':'reap','utc':now(),'pid':child.pid,'returncode':rc})
  if rc!=0: return {'ok':False,'sid':sid,'reason':f'child returncode {rc}'}
  valid_result(results/(f'roworbit-l{l}-{sid}.json'),l,a,b,source,sid)
  return {'ok':True,'sid':sid,'l':l,'start':a,'stop':b}
 except Exception as e:
  # A receipt failure after Popen is a controller failure, but never a reason
  # to orphan a child.  Wait naturally; do not signal or kill it.
  notes=[repr(e)]
  if child is not None:
   if not waited:
    try: rc=child.wait(); waited=True
    except Exception as w: notes.append('wait='+repr(w))
   for path,event in ((wait,'wait-after-receipt-failure'),(reap,'reap-after-receipt-failure')):
    if waited and not path.exists():
     try: put(path,{'event':event,'utc':now(),'pid':child.pid,'returncode':rc})
     except Exception as w: notes.append(event+'='+repr(w))
  else:
   try: put(spawn,{'event':'spawn-failed','utc':begun,'pid':None,'error':repr(e),'command':cmd})
   except Exception as w: notes.append('spawn-failed='+repr(w))
  try: put(error,{'event':'controller-receipt-or-task-error','utc':now(),'pid':None if child is None else child.pid,'waited':waited,'returncode':rc,'notes':notes})
  except Exception: pass
  return {'ok':False,'sid':sid,'reason':'; '.join(notes)}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--binary',type=Path,required=True);ap.add_argument('--source',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--workers',type=int,default=12);ap.add_argument('--chunk',type=int,default=4096);ap.add_argument('--smoke',action='store_true');ap.add_argument('--smoke-items',type=int,default=8);z=ap.parse_args()
 if z.out.exists() or not z.binary.is_file() or not z.source.is_file() or not(1<=z.workers<=12) or z.chunk<1: raise SystemExit('refusing existing output/bad args')
 if z.smoke and 'smoke' not in z.out.name.lower(): raise SystemExit('--smoke requires a distinct output directory named smoke')
 if z.smoke and z.smoke_items<1: raise SystemExit('bad smoke-items')
 source=sha(z.source); session=uuid.uuid4().hex; z.out.mkdir(); (z.out/'logs').mkdir(); (z.out/'results').mkdir()
 put(z.out/'started.json',{'schema':1,'utc':now(),'controller_pid':os.getpid(),'session':session,'N':N,'sorted_transposed_left_sizes':[7,6],'original_small_left_by_sorted_transposed_left':{'7':4,'6':5},'workers':z.workers,'chunk':z.chunk,'smoke':z.smoke,'binary':str(z.binary),'binary_sha256':sha(z.binary),'source':str(z.source),'source_sha256':source,'controller_sha256':sha(Path(__file__))})
 signal.signal(signal.SIGINT,on_signal); signal.signal(signal.SIGTERM,on_signal)
 def jobs():
  if z.smoke:
   for l in (7,6): yield l,0,min(z.smoke_items,total(l))
  else:
   for l in (7,6):
    for a in range(0,total(l),z.chunk):yield l,a,min(a+z.chunk,total(l))
 it=iter(jobs()); pending=set(); completed=0; failure=None
 with cf.ThreadPoolExecutor(max_workers=z.workers) as pool:
  while len(pending)<z.workers and not STOP:
   try: pending.add(pool.submit(task,z.binary,source,z.out,session,*next(it)))
   except StopIteration: break
  while pending:
   done,pending=cf.wait(pending,return_when=cf.FIRST_COMPLETED)
   for f in done:
    r=f.result(); completed+=1
    if not r['ok'] and failure is None: failure=r; globals()['STOP']=True; globals()['STOP_REASON']='child-or-checkpoint-failure'
    put(z.out/'progress.json',{'utc':now(),'completed_callbacks':completed,'active_callbacks':len(pending),'stop_dispatch':STOP,'stop_reason':STOP_REASON,'failure':failure})
   while len(pending)<z.workers and not STOP:
    try: pending.add(pool.submit(task,z.binary,source,z.out,session,*next(it)))
    except StopIteration: break
 status='dispatch-complete-pending-independent-audit' if not STOP and failure is None else 'stopped-after-active-reaped'
 put(z.out/'controller-finished.json',{'utc':now(),'status':status,'completed_callbacks':completed,'stop_reason':STOP_REASON,'failure':failure,'all_children_waited_and_reaped':True})
 return 0 if status.startswith('dispatch-complete') else 3
if __name__=='__main__':raise SystemExit(main())
