"""Recover interrupted batches without overwriting partial ledgers.

Completed search ledgers are replayed by the independent original verifier.
Unfinished searches are preserved under .interrupted before regeneration.
A completed batch receipt covers only this recovery's elapsed wall time.
"""
import argparse,hashlib,json,os,shutil,subprocess,sys,time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed

def atomic(p,obj):
 q=p.with_suffix(p.suffix+'.tmp');q.write_text(json.dumps(obj,indent=2)+'\n');q.replace(p)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('config',type=Path);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--workers',type=int,required=True);a=ap.parse_args()
 assert Path('/Volumes/data').is_mount();base=Path(__file__).resolve().parent
 configs=json.loads(a.config.read_text());stamp=time.strftime('%Y%m%dT%H%M%S',time.gmtime());started=time.monotonic()
 search_sha=hashlib.sha256((base/'support_search_layer.py').read_bytes()).hexdigest()
 verify_sha=hashlib.sha256((base/'verify_support_layer.py').read_bytes()).hexdigest()
 def task(c,m):
  p=a.out/f"{c['tag']}-m{m:02d}";summary=None
  if (p/'summary.json').exists():summary=json.loads((p/'summary.json').read_text())
  if summary and summary['complete']:
   assert summary['script_sha256']==search_sha
   conf=summary['config'];assert (conf['min_total'],conf['total'],conf['max_terms'],conf['min_degree'],conf['max_first_degree'])==(c['min_total'],c['total'],c['max_terms'],m,m)
   if (p/'verification.json').exists():
    v=json.loads((p/'verification.json').read_text());assert v['verified'] and v['verifier_sha256']==verify_sha
    assert v['counts']==summary['counts'] and v['ledger_sha256_uncompressed']==summary['ledger_sha256_uncompressed']
    state=json.loads((p/'job.json').read_text());assert state['status']=='verified' and state['process_exited']
    return state
   prior=p/f'recovery-prior-{stamp}';prior.mkdir()
   for name in ('job.json','verify.stdout','verify.stderr'):
    if(p/name).exists():shutil.copy2(p/name,prior/name)
   phases=['verify']
  else:
   if p.exists():
    archive=a.out/'.interrupted';archive.mkdir(exist_ok=True);p.rename(archive/f'{p.name}-{stamp}')
   p.mkdir(parents=True);phases=['search','verify']
  state={**c,'m':m,'status':'started','process_exited':False,'recovery_stamp':stamp,'reused_completed_search':phases==['verify']}
  for phase in phases:
   cmd=[sys.executable,str(base/'verify_support_layer.py'),str(p)] if phase=='verify' else [sys.executable,str(base/'support_search_layer.py'),'--total',str(c['total']),'--min-total',str(c['min_total']),'--min-degree',str(m),'--max-first-degree',str(m),'--max-terms',str(c['max_terms']),'--seconds','1e100','--out',str(p)]
   with(p/(phase+'.stdout')).open('w')as out,(p/(phase+'.stderr')).open('w')as err:
    child=subprocess.Popen(cmd,stdout=out,stderr=err,env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1','OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'})
    state.update(status='running',phase=phase,pid=child.pid,process_exited=False);atomic(p/'job.json',state)
    rc=child.wait();state.update(process_exited=True,returncode=rc)
   if rc:state['status']='unknown_failure';atomic(p/'job.json',state);return state
  state['status']='verified';atomic(p/'job.json',state);print(json.dumps(state),flush=True);return state
 tasks=[(c,m)for c in configs for m in range(1,c['total']//2+1)]
 # Prioritize incomplete verifications over fresh searches.
 tasks.sort(key=lambda cm:not (a.out/f"{cm[0]['tag']}-m{cm[1]:02d}"/'summary.json').exists())
 with ThreadPoolExecutor(max_workers=a.workers)as pool:
  jobs=[f.result()for f in as_completed([pool.submit(task,c,m)for c,m in tasks])]
 report={'workers':a.workers,'configs':configs,'jobs':jobs,'deadline_seconds':None,'elapsed_seconds':time.monotonic()-started,'elapsed_scope':'this recovery invocation only; interrupted earlier runtime not included','all_workers_exited':True,'recovery_stamp':stamp,'recovery_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 atomic(a.out/'batch.json',report)
 if any(j['status']!='verified'for j in jobs):raise SystemExit(1)
if __name__=='__main__':main()
