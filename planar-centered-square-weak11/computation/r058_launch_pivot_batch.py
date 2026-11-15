#!/usr/bin/env python3
import hashlib,json,os,subprocess,sys,time
from pathlib import Path
assert os.path.ismount('/Volumes/data')
ROOT=Path(__file__).resolve().parents[1];TMP=Path('/Volumes/data/CodexTemp/T2/planar-centered-square-weak11/r058-pivot-search');TMP.mkdir(parents=True,exist_ok=True)
COMPUTE=ROOT/'computation';OUT=COMPUTE/'results/r058-pivot-search';OUT.mkdir(parents=True,exist_ok=True)
JOBS=[(8,580801),(8,580802),(10,581001),(10,581002),(12,581201),(12,581202)];running=[];entries=[];env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',TMPDIR=str(TMP),TMP=str(TMP),TEMP=str(TMP),XDG_CACHE_HOME='/Volumes/data/CodexTemp/cache')
start=time.perf_counter()
for n,seed in JOBS:
 cmd=[sys.executable,str(COMPUTE/'r058_pivot_batch.py'),'--n',str(n),'--seed',str(seed),'--count','2000'];log=TMP/f'launch-n{n}-s{seed}.stdout';err=TMP/f'launch-n{n}-s{seed}.stderr';fo=log.open('w');fe=err.open('w');p=subprocess.Popen(cmd,stdout=fo,stderr=fe,env=env);entry={'n':n,'seed':seed,'count':2000,'pid':p.pid,'command':cmd,'stdout':str(log),'stderr':str(err)};running.append((p,fo,fe,entry));entries.append(entry)
manifest={'status':'running','jobs':entries,'driver_source_sha256':hashlib.sha256((COMPUTE/'r058_pivot_batch.py').read_bytes()).hexdigest(),'kernel_source_sha256':hashlib.sha256((COMPUTE/'r058_pivot_decomposition.cpp').read_bytes()).hexdigest()};path=OUT/'six-shard-launch.json';path.write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(manifest),flush=True)
for p,fo,fe,entry in running:
 entry['exit_code']=p.wait();fo.close();fe.close();entry['stdout_text']=Path(entry['stdout']).read_text();entry['stderr_text']=Path(entry['stderr']).read_text();result=OUT/f"n{entry['n']}-s{entry['seed']}-k3.json";entry['result']=str(result)
 if result.exists():entry['result_sha256']=hashlib.sha256(result.read_bytes()).hexdigest()
manifest['status']='complete';manifest['elapsed_seconds']=time.perf_counter()-start;path.write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(manifest,indent=2),flush=True)
assert all(e['exit_code']==0 for e in entries)
