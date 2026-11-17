#!/usr/bin/env python3
"""Unique-directory, strict-host-key remote exact histogram scan with receipts."""
import base64
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

source = Path(__file__).with_name('triangular_weight_scan.cpp').read_bytes()
remote = r'''
import base64, hashlib, json, os, platform, subprocess, sys, time, threading
from pathlib import Path
root=Path('/Users/zhengzhihao/CodexRemoteMath/a28-selector/r008/triangular-k4')
root.mkdir(parents=True,exist_ok=True)
if (root/'running.json').exists() or (root/'completed.json').exists():
    raise RuntimeError('Existing task marker: inspect, do not restart.')
source=base64.b64decode(PAYLOAD)
(root/'scan.cpp').write_bytes(source)
compile_command=['clang++','-O3','-std=c++20',str(root/'scan.cpp'),'-o',str(root/'scan')]
cc=subprocess.run(compile_command,capture_output=True,text=True)
assert cc.returncode==0,cc.stderr
pilot=subprocess.run([str(root/'scan'),'0','1048576',str(root/'pilot.json')],capture_output=True,text=True)
assert pilot.returncode==0,pilot.stderr
pilot_data=json.loads((root/'pilot.json').read_text())
assert sum(pilot_data['histogram'])==1048576 and pilot_data['histogram'][0]==1
jobs=[];workers=14;total=1<<31;began=time.monotonic()
for i in range(workers):
    path=root/('shard-%02d.json'%i)
    start=total*i//workers;end=total*(i+1)//workers
    command=['nice','-n','10',str(root/'scan'),str(start),str(end-start),str(path)]
    proc=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    jobs.append((proc,command,path))
pids=[job[0].pid for job in jobs]
(root/'running.json').write_text(json.dumps({'pids':pids,'source_sha256':hashlib.sha256(source).hexdigest()}))
print(json.dumps({'event':'started','pids':pids}),flush=True)
cpu=[];top={}
def sample_top():
    p=subprocess.run(['top','-l','2','-s','1','-n','0'],capture_output=True,text=True)
    top.update(returncode=p.returncode,lines=[s for s in p.stdout.splitlines() if 'CPU usage' in s])
th=threading.Thread(target=sample_top);th.start()
while any(p.poll() is None for p,_,_ in jobs):
    p=subprocess.run(['ps','-p',','.join(map(str,pids)),'-o','pid=,pcpu=,state='],capture_output=True,text=True)
    cpu.append({'seconds':time.monotonic()-began,'percent_sum':sum(float(s.split()[1]) for s in p.stdout.splitlines() if len(s.split())>=2)})
    time.sleep(0.2)
th.join();shards=[]
for p,command,path in jobs:
    stdout,stderr=p.communicate()
    assert p.returncode==0,(p.pid,p.returncode,stderr)
    shards.append({'pid':p.pid,'command':command,'exit_code':p.returncode,'stdout':stdout,'stderr':stderr,'result':json.loads(path.read_text())})
assert sum(s['result']['count'] for s in shards)==total
result={'status':'passed','host':platform.node(),'platform':platform.platform(),'python':sys.version,
        'logical_cpus':os.cpu_count(),'compiler':subprocess.check_output(['clang++','--version'],text=True),
        'compile_command':compile_command,'source_sha256':hashlib.sha256(source).hexdigest(),
        'pilot':pilot_data,'workers':workers,'shards':shards,'cpu_samples':cpu,'whole_machine_top':top,
        'elapsed_seconds':time.monotonic()-began}
(root/'completed.json').write_text(json.dumps(result,indent=2)+'\n')
(root/'running.json').unlink()
print(json.dumps({'event':'completed','result':result}),flush=True)
'''.replace('PAYLOAD', repr(base64.b64encode(source).decode()))
start=time.monotonic()
p=subprocess.run(['ssh','-o','BatchMode=yes','-o','StrictHostKeyChecking=yes','-o','ConnectTimeout=15',
                  'macbook-math','/usr/bin/python3 -'],input=remote,capture_output=True,text=True)
receipt={'ssh_exit_code':p.returncode,'elapsed_seconds':time.monotonic()-start,
         'stdout':p.stdout,'stderr':p.stderr,'runner_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
for line in p.stdout.splitlines():
    if line.startswith('{'):
        event=json.loads(line)
        if event['event']=='completed':receipt['result']=event['result']
Path(sys.argv[1]).write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'ssh_exit_code':p.returncode,'complete':'result' in receipt,'elapsed_seconds':receipt['elapsed_seconds']}))
raise SystemExit(p.returncode)
