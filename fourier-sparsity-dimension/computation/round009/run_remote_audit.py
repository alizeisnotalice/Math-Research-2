#!/usr/bin/env python3
"""Run one identified exact certificate replay on the authorized MacBook."""
import argparse,base64,hashlib,json,subprocess,time
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args()
base=Path(__file__).parent
source=base.joinpath('recursive_normal_rank_certificate.py').read_bytes()
input_data=base.parent.joinpath('round008/remote-packing-execution.json').read_bytes()
remote=r'''
import base64,hashlib,json,os,platform,subprocess,sys,time
from pathlib import Path
root=Path('/Users/zhengzhihao/CodexRemoteMath/a28-selector/r009/recursive-audit')
root.mkdir(parents=True,exist_ok=False)
source=base64.b64decode(SOURCE);input_data=base64.b64decode(INPUT)
(root/'verify.py').write_bytes(source);(root/'packing-input.json').write_bytes(input_data)
command=['nice','-n','10','/usr/bin/python3',str(root/'verify.py'),'--out',str(root/'result.json'),'--remote-packing',str(root/'packing-input.json')]
p=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
(root/'running.json').write_text(json.dumps({'pid':p.pid,'command':command}))
beg=time.monotonic();stdout,stderr=p.communicate()
receipt={'host':platform.node(),'platform':platform.platform(),'python':sys.version,'logical_cpus':os.cpu_count(),
 'pid':p.pid,'command':command,'exit_code':p.returncode,'stdout':stdout,'stderr':stderr,
 'source_sha256':hashlib.sha256(source).hexdigest(),'input_sha256':hashlib.sha256(input_data).hexdigest(),
 'elapsed_seconds':time.monotonic()-beg}
if p.returncode==0:receipt['result']=json.loads((root/'result.json').read_text())
(root/'completed.json').write_text(json.dumps(receipt,indent=2)+'\n');(root/'running.json').unlink()
print(json.dumps(receipt));sys.exit(p.returncode)
'''.replace('SOURCE',repr(base64.b64encode(source).decode())).replace('INPUT',repr(base64.b64encode(input_data).decode()))
beg=time.monotonic()
p=subprocess.run(['ssh','-o','BatchMode=yes','-o','StrictHostKeyChecking=yes','-o','ConnectTimeout=10','macbook-math','/usr/bin/python3 -'],input=remote,capture_output=True,text=True)
receipt={'ssh_exit_code':p.returncode,'elapsed_seconds':time.monotonic()-beg,'stderr':p.stderr,'runner_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
if p.stdout.strip():receipt['remote']=json.loads(p.stdout)
a.out.write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'ssh_exit_code':p.returncode,'elapsed_seconds':receipt['elapsed_seconds'],'remote_status':receipt.get('remote',{}).get('exit_code')}))
raise SystemExit(p.returncode)
