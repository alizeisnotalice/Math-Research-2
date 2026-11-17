#!/Users/zzh/.venvs/main/bin/python
"""Controller only; execution requires explicit --execute after parent signal."""
import argparse,concurrent.futures,hashlib,json,subprocess,sys,time
from pathlib import Path
H=Path(__file__).resolve().parent; SOURCE=H/'s1_lowprefix_screen.py'
def tn():
 from math import isqrt
 for T in range(-16,17):
  for D in range(1,(17-abs(T))**2):
   if (T*T-D)%4==0 and isqrt(D)**2!=D:yield T,(T*T-D)//4,D
def work(ix):
 t=list(tn())[ix]; cmd=[sys.executable,str(SOURCE),'--tn-start',str(ix),'--tn-count','1']
 p=subprocess.run(cmd,capture_output=True,text=True,timeout=4*3600)
 raw=H/f's1_shard_{ix}_1.json'; rec={'index':ix,'TNDelta':t,'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'returncode':p.returncode,'stdout':p.stdout,'stderr':p.stderr}
 if p.returncode==0 and raw.exists(): rec['result']=json.loads(raw.read_text());rec['status']='completed'
 else: rec['status']='exception'
 out=H/'s1_tn_records';out.mkdir(exist_ok=True); (out/f'{ix:04d}.json').write_text(json.dumps(rec,indent=2)+'\n');return rec['status']
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--execute',action='store_true');ap.add_argument('--workers',type=int,default=4);z=ap.parse_args(); rows=list(tn()); manifest={'status':'PREPARED_NOT_RUN','TN_total':len(rows),'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'workers':z.workers,'execute':z.execute,'submitted':list(range(len(rows)))}
 if z.execute:
  with concurrent.futures.ThreadPoolExecutor(max_workers=z.workers) as ex: manifest['completed']=[f.result() for f in concurrent.futures.as_completed([ex.submit(work,i) for i in range(len(rows))])]
  manifest['status']='FINISHED'
 (H/'s1_controller_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps({'status':manifest['status'],'TN_total':len(rows)}))
if __name__=='__main__':main()
