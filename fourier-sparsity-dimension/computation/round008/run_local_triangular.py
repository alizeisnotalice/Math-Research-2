#!/usr/bin/env python3
"""Bounded local fallback for the complete triangular weight scan."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

ap=argparse.ArgumentParser();ap.add_argument('--work',type=Path,required=True)
ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
assert Path('/Volumes/data').is_mount()
a.work.mkdir(parents=True,exist_ok=True)
if (a.work/'running.json').exists() or (a.work/'completed.json').exists():
    raise RuntimeError('Existing task: inspect before replay.')
here=Path(__file__).resolve().parent;source=here/'triangular_weight_scan.cpp'
binary=a.work/'scan';env=os.environ.copy();env['CLANG_MODULE_CACHE_PATH']='/Volumes/data/CodexTemp/cache/clang'
cmd=['clang++','-O3','-std=c++20',str(source),'-o',str(binary)]
subprocess.run(cmd,check=True,env=env,capture_output=True,text=True)
workers=min(8,os.cpu_count() or 1);total=1<<31;jobs=[];began=time.monotonic()
for i in range(workers):
    start=total*i//workers;end=total*(i+1)//workers;path=a.work/f'shard-{i:02d}.json'
    command=['nice','-n','10',str(binary),str(start),str(end-start),str(path)]
    proc=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    jobs.append((proc,command,path))
(a.work/'running.json').write_text(json.dumps({'pids':[p.pid for p,_,_ in jobs]}))
shards=[]
for p,command,path in jobs:
    stdout,stderr=p.communicate();assert p.returncode==0,(p.returncode,stderr)
    shards.append({'pid':p.pid,'command':command,'exit_code':p.returncode,'stdout':stdout,'stderr':stderr,
                   'result':json.loads(path.read_text())})
result={'status':'passed','host':platform.node(),'platform':platform.platform(),'python':sys.version,
        'logical_cpus':os.cpu_count(),'workers':workers,'compile_command':cmd,
        'compiler':subprocess.check_output(['clang++','--version'],text=True),
        'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
        'binary_sha256':hashlib.sha256(binary.read_bytes()).hexdigest(),
        'shards':shards,'elapsed_seconds':time.monotonic()-began}
receipt={'execution_exit_code':0,'execution':'local CPU fallback; remote connection unavailable',
         'result':result,'runner_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
a.out.write_text(json.dumps(receipt,indent=2)+'\n')
(a.work/'completed.json').write_text(json.dumps({'receipt':str(a.out),'status':'passed'}))
(a.work/'running.json').unlink()
print(json.dumps({'status':'passed','workers':workers,'tested':sum(s['result']['count'] for s in shards),'seconds':result['elapsed_seconds']}))
