import subprocess,json,os,sys,time,hashlib
from pathlib import Path
r=Path('/Users/zhengzhihao/CodexRemoteMath/a30-semantic51-52')
(r/'supervisor-pid.json').write_text(json.dumps({'pid':os.getpid(),'start':time.time()}))
for total in (51,52):
 cmd=[sys.executable,str(r/'code/runner.py'),'--root',str(r),'--total',str(total),'--workers','13','--seconds','14400']
 with (r/f'compute{total}.log').open('a') as log:
  p=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1','OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'})
  (r/f'launch{total}.json').write_text(json.dumps({'pid':p.pid,'supervisor':os.getpid(),'command':cmd,'time':time.time(),'runner_sha256':hashlib.sha256((r/'code/runner.py').read_bytes()).hexdigest()}))
  rc=p.wait()
 summary=json.loads((r/f'progress{total}.json').read_text()) if (r/f'progress{total}.json').exists() else {}
 (r/f'exit{total}.json').write_text(json.dumps({'returncode':rc,'pid':p.pid,'child_waited':True,'time':time.time(),'status':summary.get('status'),'verified':summary.get('verified'),'remaining':summary.get('remaining')}))
 if rc or summary.get('status')!='complete':raise SystemExit(rc or 3)
