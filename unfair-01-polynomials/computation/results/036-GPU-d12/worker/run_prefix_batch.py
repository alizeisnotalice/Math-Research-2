"""Run one owned, checkpointing MPS batch and wait/reap before exit."""
import argparse,hashlib,json,os,signal,subprocess,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
WORKER=HERE/'gpu_factor_search_prefix.py'
def atom(p,v):
 t=p.with_suffix('.tmp');t.write_text(json.dumps(v,indent=2)+'\n');t.replace(p)
def stop(child):
 for sig in [signal.SIGINT,signal.SIGTERM,signal.SIGKILL]:
  if child.poll() is not None:return
  os.killpg(child.pid,sig)
  try:child.wait(timeout=120);return
  except subprocess.TimeoutExpired:pass
 child.wait()
def main():
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True)
 for key,default in [('m',32),('n',72),('d',7),('starts',16384),('steps',12000),('seed',2026111605),('seconds',14400),('restart',2000)]:p.add_argument('--'+key,type=int,default=default)
 p.add_argument('--margin',type=float,default=.001);a=p.parse_args()
 assert Path('/Volumes/data').is_mount() and a.seconds>0
 assert not a.out.exists() or not any(a.out.iterdir());a.out.mkdir(parents=True,exist_ok=True)
 assert hashlib.sha256(WORKER.read_bytes()).hexdigest()=='250ecf20c4ac60a67385c4e3e64cd856639aee543e2cdbc11ae96e7630ca2cad'
 cmd=[sys.executable,str(WORKER),'--m',str(a.m),'--n',str(a.n),'--fractional-index',str(a.d),'--starts',str(a.starts),'--steps',str(a.steps),'--seed',str(a.seed),'--margin',str(a.margin),'--restart-every',str(a.restart),'--out',str(a.out/'run')]
 log=(a.out/'worker.log').open('w');child=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,start_new_session=True);start=time.monotonic()
 record={'status':'running','parent_pid':os.getpid(),'child_pid':child.pid,'command':cmd,'seconds_limit':a.seconds,'worker_sha256':hashlib.sha256(WORKER.read_bytes()).hexdigest(),'driver_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'All filtered legal low prefixes for this d are assigned among starts; numerical trajectories only, not an exclusion.'};atom(a.out/'queue.json',record);interrupted=False
 try:
  while child.poll() is None and time.monotonic()-start<a.seconds:time.sleep(5)
  if child.poll() is None:record['status']='partial_time_limit';stop(child)
  else:record['status']='complete' if child.returncode==0 else 'failed'
 except KeyboardInterrupt:
  interrupted=True;record['status']='partial_user_or_driver_interrupt';stop(child)
 finally:
  if child.poll() is None:stop(child)
  child.wait();log.close();record.update(returncode=child.returncode,elapsed_seconds=time.monotonic()-start,child_reaped=True,release_exists=(a.out/'run/release.json').is_file());atom(a.out/'queue.json',record)
 print(json.dumps(record))
 if interrupted:raise SystemExit(130)
 if record['status']!='complete':raise SystemExit(1)
if __name__=='__main__':main()
