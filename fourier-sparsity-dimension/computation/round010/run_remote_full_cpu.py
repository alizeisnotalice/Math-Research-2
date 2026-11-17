#!/usr/bin/env python3
"""Bounded 15-process exact Walsh audits on the authorized MacBook."""
import argparse,base64,hashlib,json,subprocess,time
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args()
source=Path(__file__).with_name('cddl_fact22_lift_verify.cpp').read_bytes()
remote=r'''
import base64,hashlib,json,os,platform,subprocess,sys,threading,time
from pathlib import Path
root=Path('/Users/zhengzhihao/CodexRemoteMath/a28-selector/r010/lift-15-full')
root.mkdir(parents=True,exist_ok=False)
source=base64.b64decode(PAYLOAD);(root/'verify.cpp').write_bytes(source)
cache=root.parent.parent/'cache/clang';cache.mkdir(parents=True,exist_ok=True)
env=os.environ.copy();env.update(TMPDIR=str(root),TMP=str(root),TEMP=str(root))
compile_command=['clang++','-O3','-std=c++20','-fmodules-cache-path='+str(cache),str(root/'verify.cpp'),'-o',str(root/'verify')]
compiled=subprocess.run(compile_command,capture_output=True,text=True,env=env)
assert compiled.returncode==0,compiled.stderr
pilot_command=[str(root/'verify'),'--max-l','1','--samples','1','--seed','202611172000','--out',str(root/'pilot.json')]
pilot=subprocess.run(pilot_command,capture_output=True,text=True);assert pilot.returncode==0,pilot.stderr
jobs=[];began=time.monotonic()
for i in range(15):
 path=root/('job-%02d.json'%i)
 cmd=[str(root/'verify'),'--max-l','4','--samples','1024','--seed',str(202611172000+i),'--out',str(path)]
 proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
 jobs.append((proc,cmd,path))
pids=[p.pid for p,_,_ in jobs]
(root/'running.json').write_text(json.dumps({'pids':pids,'source_sha256':hashlib.sha256(source).hexdigest()}))
print(json.dumps({'event':'started','workers':15,'pids':pids}),flush=True)
cpu=[];top={}
def sample_top():
 p=subprocess.run(['top','-l','2','-s','1','-n','0'],capture_output=True,text=True)
 top.setdefault('samples',[]).append({'seconds':time.monotonic()-began,'returncode':p.returncode,'lines':[x for x in p.stdout.splitlines() if 'CPU usage' in x]})
def monitor_top():
 while any(p.poll() is None for p,_,_ in jobs):
  sample_top();time.sleep(5)
th=threading.Thread(target=monitor_top);th.start();last=-10
while any(p.poll() is None for p,_,_ in jobs):
 p=subprocess.run(['ps','-p',','.join(map(str,pids)),'-o','pid=,pcpu=,state='],capture_output=True,text=True)
 sample={'seconds':time.monotonic()-began,'processes':p.stdout.splitlines(),'percent_sum':sum(float(x.split()[1]) for x in p.stdout.splitlines() if len(x.split())>=2),'active_workers':sum(p.poll() is None for p,_,_ in jobs)}
 cpu.append(sample)
 if sample['seconds']-last>=10:print(json.dumps({'event':'cpu',**sample}),flush=True);last=sample['seconds']
 time.sleep(.5)
th.join();records=[]
for p,cmd,path in jobs:
 stdout,stderr=p.communicate();assert p.returncode==0,(p.pid,stderr)
 result=json.loads(path.read_text());assert result['status']=='passed'
 samples=result.pop('samples');assert len(samples)==5120
 assert all(all(s[k] for k in ('pi_bijective','all_walsh_abs_q','direct_walsh_matches_formula','one_bit_flip_breaks_flatness')) for s in samples)
 result['sample_count']=len(samples);result['direct_coefficient_count']=sum(1<<s['variables'] for s in samples)
 result['raw_result_sha256']=hashlib.sha256(path.read_bytes()).hexdigest();result['raw_result_path']=str(path)
 result['all_sample_predicates_verified']=True
 records.append({'pid':p.pid,'command':cmd,'exit_code':p.returncode,'stdout':stdout,'stderr':stderr,'result':result})
result={'status':'passed','host':platform.node(),'platform':platform.platform(),'python':sys.version,'logical_cpus':os.cpu_count(),
 'compiler':subprocess.check_output(['clang++','--version'],text=True),'compile_command':compile_command,
 'source_sha256':hashlib.sha256(source).hexdigest(),'binary_sha256':hashlib.sha256((root/'verify').read_bytes()).hexdigest(),
 'pilot':{'command':pilot_command,'exit_code':pilot.returncode,'result':json.loads((root/'pilot.json').read_text())},
 'jobs':records,'cpu_samples':cpu,'whole_machine_top':top,'elapsed_seconds':time.monotonic()-began}
(root/'completed.json').write_text(json.dumps(result,indent=2)+'\n');(root/'running.json').unlink()
print(json.dumps({'event':'completed','result':result}),flush=True)
'''.replace('PAYLOAD',repr(base64.b64encode(source).decode()))
began=time.monotonic();p=subprocess.Popen(['ssh','-o','BatchMode=yes','-o','StrictHostKeyChecking=yes','-o','ConnectTimeout=10','macbook-math','/usr/bin/python3 -'],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
p.stdin.write(remote);p.stdin.close();result=None
for line in p.stdout:
 if line.startswith('{'):
  event=json.loads(line)
  if event['event']=='completed':result=event['result']
  else:print(line.strip(),flush=True)
code=p.wait();stderr=p.stderr.read()
receipt={'ssh_exit_code':code,'elapsed_seconds':time.monotonic()-began,'stderr':stderr,'runner_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
if result is not None:receipt['result']=result
a.out.write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'ssh_exit_code':code,'complete':result is not None,'elapsed_seconds':receipt['elapsed_seconds']}));raise SystemExit(code)
