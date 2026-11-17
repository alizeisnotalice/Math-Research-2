#!/usr/bin/env python3
"""Strict v2 GPU52 bridge; v1 records remain a separately audited layer."""
from __future__ import annotations
import argparse,base64,hashlib,json,os,subprocess,time
from pathlib import Path

SSH=["-o","HostName=fe80::8c00:aeff:feee:2fe%%anpi1","-o","HostKeyAlias=169.254.63.75","-o","ControlMaster=auto","-o","ControlPersist=600","-o","ControlPath=/Volumes/data/CodexTemp/T2/unfair-01-polynomials/macbook-anpi1-20261117.sock","-o","BatchMode=yes","-o","StrictHostKeyChecking=yes","-o","ConnectTimeout=10"]
CMD="nice -n 10 env PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /Users/zhengzhihao/.venvs/math/bin/python -u -"
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def sync(f):f.flush();os.fsync(f.fileno())
def syncd(p):
 fd=os.open(p,os.O_RDONLY)
 try:os.fsync(fd)
 finally:os.close(fd)
def reap(p):
 if p.poll() is None:
  p.terminate()
  try:return p.wait(10)
  except subprocess.TimeoutExpired:p.kill();return p.wait(10)
 return p.returncode
def validate(data):
 assert [t.get('id') for t in data['tasks']]==list(range(24));starts={0:67108864,1:67108864,2:23166976}
 for t in data['tasks']:
  assert t['resume_start_mask']==starts.get(t['id'],0) and 0<=t['resume_start_mask']<=t['masks']
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--search',type=Path,required=True);ap.add_argument('--remote',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);mode=ap.add_mutually_exclusive_group();mode.add_argument('--run',action='store_true');mode.add_argument('--smoke',action='store_true');a=ap.parse_args()
 if not Path('/Volumes/data').is_mount() or a.out.exists():raise SystemExit('data unavailable or output exists')
 data=json.loads(a.input.read_text());validate(data);payload='INPUT_B64='+repr(base64.b64encode(a.input.read_bytes()).decode())+'\nSEARCH_B64='+repr(base64.b64encode(a.search.read_bytes()).decode())+'\n'+a.remote.read_text();a.out.mkdir(parents=True);syncd(a.out.parent)
 with (a.out/'payload.py').open('x') as f:f.write(payload);sync(f)
 syncd(a.out); manifest={'format':'a30-gpu52-resume-v2','mode':'smoke' if a.smoke else 'run' if a.run else 'dry_run','scope':'Original received GPU52 v1 prefixes used the correct formula and remain a separate audited layer. This v2 continuation begins only at v1 received boundaries; no global completion claim.','input_sha256':sha(a.input),'search_sha256':sha(a.search),'remote_sha256':sha(a.remote),'payload_sha256':hashlib.sha256(payload.encode()).hexdigest(),'task_starts':{str(t['id']):t['resume_start_mask'] for t in data['tasks']},'workers':1,'seconds_limit':14400,'ssh_command':['ssh',*SSH,'macbook-math',CMD]}
 with (a.out/'launch-manifest.json').open('x') as f:json.dump(manifest,f,indent=2,sort_keys=True);f.write('\n');sync(f)
 if not(a.run or a.smoke):print(json.dumps(manifest,sort_keys=True));return
 command=['ssh',*SSH,'macbook-math','env','GPU52_SMOKE='+('1' if a.smoke else '0'),CMD];proc=None;failure=None;rc=None;lines=batches=0
 with (a.out/'stderr.txt').open('x') as err,(a.out/'stdout.jsonl').open('x',buffering=1) as out:
  try:
   proc=subprocess.Popen(command,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=err,text=True,bufsize=1);assert proc.stdin and proc.stdout;proc.stdin.write(payload);proc.stdin.close()
   for line in proc.stdout:
    out.write(line);sync(out);lines+=1
    if line.startswith('{"event":"batch"'):batches+=1
   rc=proc.wait()
  except BaseException as e:
   failure={'type':type(e).__name__,'message':str(e)}
   if proc is not None:rc=reap(proc)
  finally:
   if proc is not None and proc.poll() is None:rc=reap(proc)
   sync(err)
 receipt={'returncode':rc,'failure':failure,'ssh_pid':None if proc is None else proc.pid,'lines_received':lines,'batches_received':batches,'all_owned_processes_reaped':proc is None or proc.poll() is not None,'finished_unix':time.time(),**manifest}
 with (a.out/'exit.json').open('x') as f:json.dump(receipt,f,indent=2,sort_keys=True);f.write('\n');sync(f)
 raise SystemExit(0 if rc==0 and failure is None else 1)
if __name__=='__main__':main()
