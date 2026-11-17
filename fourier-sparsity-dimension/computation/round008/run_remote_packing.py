#!/usr/bin/env python3
"""Fourteen independent, finite, calibrated MacBook packing-search jobs."""
import base64
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

source=Path(__file__).with_name('packing_search.cpp').read_bytes()
remote=r'''
import base64,hashlib,json,os,platform,subprocess,sys,time,threading
from pathlib import Path
root=Path('/Users/zhengzhihao/CodexRemoteMath/a28-selector/r008/packing-14')
root.mkdir(parents=True,exist_ok=True)
if (root/'running.json').exists() or (root/'completed.json').exists():
    raise RuntimeError('Existing task record; inspect without restarting.')
source=base64.b64decode(PAYLOAD);(root/'search.cpp').write_bytes(source)
compile_command=['clang++','-O3','-std=c++20',str(root/'search.cpp'),'-o',str(root/'search')]
subprocess.run(compile_command,check=True,capture_output=True,text=True)
pilot_records=[];iteration_budget={}
for n in (8,9):
    path=root/('pilot-n%d.json'%n)
    command=[str(root/'search'),'--n',str(n),'--seed',str(202611170800+n),
             '--iterations','1000000','--target','8','--out',str(path)]
    p=subprocess.run(command,capture_output=True,text=True);assert p.returncode==0,p.stderr
    data=json.loads(path.read_text());assert data['status']=='passed'
    if data['counterexample_found']:
        (root/'pilot-counterexample.json').write_text(json.dumps(data,indent=2))
        print(json.dumps({'event':'pilot_counterexample','result':data}),flush=True)
        raise SystemExit(4)
    # Target about 90 single-process seconds, bounded to at most 100M iterations.
    iterations=max(1000000,min(100000000,int(90*1000000/max(data['elapsed_seconds'],0.01))))
    iteration_budget[n]=iterations
    pilot_records.append({'command':command,'exit_code':p.returncode,'result':data})
jobs=[];began=time.monotonic()
for i in range(14):
    n=8+(i%2);seed=2026111708000+i;path=root/('job-%02d.json'%i)
    cmd=['nice','-n','10',str(root/'search'),'--n',str(n),'--seed',str(seed),
         '--iterations',str(iteration_budget[n]),'--target','8','--out',str(path)]
    p=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    jobs.append((p,cmd,path))
pids=[p.pid for p,_,_ in jobs]
(root/'running.json').write_text(json.dumps({'pids':pids,'started':time.time(),'source_sha256':hashlib.sha256(source).hexdigest()}))
print(json.dumps({'event':'started','workers':len(jobs),'pids':pids,'iterations_by_n':iteration_budget}),flush=True)
cpu=[];top={}
def sample_top():
    p=subprocess.run(['top','-l','2','-s','1','-n','0'],capture_output=True,text=True)
    top.update(returncode=p.returncode,lines=[line for line in p.stdout.splitlines() if 'CPU usage' in line])
th=threading.Thread(target=sample_top);th.start();last=-10
while any(p.poll() is None for p,_,_ in jobs):
    p=subprocess.run(['ps','-p',','.join(map(str,pids)),'-o','pid=,pcpu=,state='],capture_output=True,text=True)
    sample={'seconds':time.monotonic()-began,'processes':p.stdout.splitlines(),
            'percent_sum':sum(float(line.split()[1]) for line in p.stdout.splitlines() if len(line.split())>=2),
            'active_workers':sum(j[0].poll() is None for j in jobs)}
    cpu.append(sample)
    if sample['seconds']-last>=5:
        print(json.dumps({'event':'cpu','seconds':sample['seconds'],'percent_sum':sample['percent_sum'],
                          'active_workers':sample['active_workers']}),flush=True);last=sample['seconds']
    time.sleep(0.5)
th.join();records=[]
for p,cmd,path in jobs:
    stdout,stderr=p.communicate();assert p.returncode==0,(p.pid,p.returncode,stderr)
    result=json.loads(path.read_text());assert result['status']=='passed'
    records.append({'pid':p.pid,'command':cmd,'exit_code':p.returncode,'stdout':stdout,'stderr':stderr,'result':result})
result={'status':'passed','host':platform.node(),'platform':platform.platform(),'python':sys.version,
        'logical_cpus':os.cpu_count(),'compiler':subprocess.check_output(['clang++','--version'],text=True),
        'compile_command':compile_command,'source_sha256':hashlib.sha256(source).hexdigest(),
        'pilot_records':pilot_records,'workers':14,'jobs':records,'cpu_samples':cpu,'whole_machine_top':top,
        'elapsed_seconds':time.monotonic()-began}
(root/'completed.json').write_text(json.dumps(result,indent=2)+'\n');(root/'running.json').unlink()
print(json.dumps({'event':'completed','result':result}),flush=True)
'''.replace('PAYLOAD',repr(base64.b64encode(source).decode()))
began=time.monotonic()
p=subprocess.Popen(['ssh','-o','BatchMode=yes','-o','StrictHostKeyChecking=yes','-o','ConnectTimeout=15',
                    'macbook-math','/usr/bin/python3 -'],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
p.stdin.write(remote);p.stdin.close();lines=[];result=None
for line in p.stdout:
    lines.append(line)
    if line.startswith('{'):
        event=json.loads(line)
        if event['event']=='completed':result=event['result']
        else:print(line.strip(),flush=True)
code=p.wait();stderr=p.stderr.read()
receipt={'ssh_exit_code':code,'elapsed_seconds':time.monotonic()-began,'stdout':''.join(lines),'stderr':stderr,
         'runner_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
if result is not None:receipt['result']=result
Path(sys.argv[1]).write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'event':'receipt_saved','ssh_exit_code':code,'complete':result is not None,
                  'elapsed_seconds':receipt['elapsed_seconds']}),flush=True)
raise SystemExit(code)
