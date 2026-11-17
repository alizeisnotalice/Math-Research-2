#!/usr/bin/env python3
"""Run disjoint exact RM3/W shards on the authorized MacBook, saving receipts."""
import base64
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
payload = {n: base64.b64encode((HERE/n).read_bytes()).decode() for n in ['rm3_multineighbor.cpp','neighbor_input.hpp']}
remote = r'''
import base64, hashlib, json, os, platform, subprocess, sys, time, threading
from pathlib import Path
payload = PAYLOAD
root=Path('/Users/zhengzhihao/CodexRemoteMath/a28-selector/r007/neighbor-keep-4567')
root.mkdir(parents=True,exist_ok=True)
# This unique directory must have no previously active or completed full run.
if (root/'full-receipt.json').exists() or (root/'running.json').exists():
    raise RuntimeError('existing run record; inspect before resuming or replacing')
for name,data in payload.items(): (root/name).write_bytes(base64.b64decode(data))
hashes={n:hashlib.sha256((root/n).read_bytes()).hexdigest() for n in payload}
cc=subprocess.run(['clang++','-O3','-std=c++20',str(root/'rm3_multineighbor.cpp'),'-o',str(root/'search')],capture_output=True,text=True)
assert cc.returncode==0,cc.stderr
compiler=subprocess.check_output(['clang++','--version'],text=True).splitlines()[0]
pilot=subprocess.run(['nice','-n','10',str(root/'search'),'0','1048576',str(root/'pilot')],capture_output=True,text=True)
assert pilot.returncode==0,pilot.stderr
pilot_summary=json.loads((root/'pilot.summary.json').read_text())
pilot_hits=[json.loads(l) for l in (root/'pilot.found.jsonl').read_text().splitlines()]
assert any(h['index']==0 for h in pilot_hits)
print(json.dumps({'event':'pilot_complete','summary':pilot_summary}),flush=True)
workers=14; total=1<<31
jobs=[]; began=time.monotonic()
for i in range(workers):
    start=total*i//workers; end=total*(i+1)//workers
    prefix=root/('shard-%02d'%i)
    out=open(str(prefix)+'.stdout','w');err=open(str(prefix)+'.stderr','w')
    cmd=['nice','-n','10',str(root/'search'),str(start),str(end-start),str(prefix)]
    proc=subprocess.Popen(cmd,stdout=out,stderr=err)
    jobs.append((proc,prefix,out,err,cmd))
pids=[j[0].pid for j in jobs]
(root/'running.json').write_text(json.dumps({'pids':pids,'hashes':hashes,'started':time.time()}))
print(json.dumps({'event':'workers_started','pids':pids,'workers':workers}),flush=True)
cpu_samples=[]
top_holder={}
def sample_top():
    top_holder['returncode']=None
    pr=subprocess.run(['top','-l','2','-s','1','-n','0'],capture_output=True,text=True)
    top_holder.update(returncode=pr.returncode,lines=[l for l in pr.stdout.splitlines() if 'CPU usage' in l])
th=threading.Thread(target=sample_top);th.start()
while any(j[0].poll() is None for j in jobs):
    pr=subprocess.run(['ps','-p',','.join(map(str,pids)),'-o','pid=,pcpu=,state='],capture_output=True,text=True)
    lines=pr.stdout.splitlines()
    cpu_samples.append({'seconds':time.monotonic()-began,'processes':lines,'percent_sum':sum(float(l.split()[1]) for l in lines if len(l.split())>=2)})
    time.sleep(0.2)
th.join()
records=[];hits=[]
for proc,prefix,out,err,cmd in jobs:
    code=proc.wait();out.close();err.close()
    assert code==0,(proc.pid,code,(Path(str(prefix)+'.stderr').read_text()))
    summary=json.loads(Path(str(prefix)+'.summary.json').read_text())
    found=[json.loads(l) for l in Path(str(prefix)+'.found.jsonl').read_text().splitlines()]
    assert summary['hits']==len(found)
    records.append({'pid':proc.pid,'exit_code':code,'command':cmd,'summary':summary})
    hits.extend(found)
assert sum(x['summary']['tested'] for x in records)==total
result={'status':'passed','host':platform.node(),'python':sys.version,'platform':platform.platform(),
        'logical_cpus':os.cpu_count(),'compiler':compiler,'source_hashes':hashes,'workers':workers,
        'pilot':pilot_summary,'shards':records,'found':hits,'total_tested':total,
        'cpu_samples':cpu_samples,'whole_machine_top':top_holder,'elapsed_seconds':time.monotonic()-began}
(root/'full-receipt.json').write_text(json.dumps(result,indent=2)+'\n')
(root/'running.json').unlink()
print(json.dumps({'event':'full_result','result':result}),flush=True)
'''.replace('PAYLOAD',repr(payload))
out_path=Path(sys.argv[1]);began=time.monotonic()
pr=subprocess.run(['ssh','-o','BatchMode=yes','-o','StrictHostKeyChecking=yes','-o','ConnectTimeout=15','macbook-math','/usr/bin/python3 -'],input=remote,capture_output=True,text=True)
receipt={'ssh_exit_code':pr.returncode,'elapsed_seconds':time.monotonic()-began,'stdout':pr.stdout,'stderr':pr.stderr,'runner_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
events=[json.loads(l) for l in pr.stdout.splitlines() if l.startswith('{')]
for e in events:
    if e['event']=='full_result':receipt['result']=e['result']
out_path.write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'ssh_exit_code':pr.returncode,'events':[e['event'] for e in events],'receipt':str(out_path),'elapsed_seconds':receipt['elapsed_seconds']}))
if 'result' in receipt:
    r=receipt['result'];print(json.dumps({'total_tested':r['total_tested'],'hits':len(r['found']),'workers':r['workers'],'max_A28_process_cpu_percent':max([s['percent_sum'] for s in r['cpu_samples']]+[0]),'whole_machine_top':r['whole_machine_top']}))
raise SystemExit(pr.returncode)
