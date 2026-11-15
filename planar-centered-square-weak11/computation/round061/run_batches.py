import json, os, subprocess, sys, time
from pathlib import Path
assert os.path.ismount('/Volumes/data')
root=Path(__file__).resolve().parent
env=dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
pending=list(range(12)); active={}; done=[]
while pending or active:
    while pending and len(active)<10:
        i=pending.pop(0); log=(root/f'batch-{i}.log').open('w')
        cmd=[sys.executable,'-B',str(root/'finite_source_cut_certificate.py'),'--batch',str(i),'--bits','24','--input',str(root/'r061-bounded-inputs.json'),'--out',str(root/f'batch-{i}.json')]
        p=subprocess.Popen(cmd,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT)
        active[p.pid]=(p,log,i,time.time(),cmd); print(f'START batch={i} pid={p.pid}',flush=True)
    for pid,(p,log,i,start,cmd) in list(active.items()):
        rc=p.poll()
        if rc is None: continue
        log.close(); done.append(dict(batch=i,pid=pid,exit_status=rc,wall_seconds=time.time()-start,argv=cmd)); print(f'END batch={i} pid={pid} exit={rc}',flush=True); del active[pid]
    time.sleep(.1)
(root/'run-receipt.json').write_text(json.dumps(done,indent=2)+'\n')
if any(x['exit_status'] for x in done): raise SystemExit(1)
