import base64,hashlib,json,os,subprocess,time
from pathlib import Path
assert os.path.ismount('/Volumes/data')
r=Path(__file__).resolve().parent
payload='INPUT_B64='+repr(base64.b64encode((r/'inputs.json').read_bytes()).decode())+'\nSEARCH_B64='+repr(base64.b64encode((r/'search.wl').read_bytes()).decode())+'\n'+(r/'remote_gpu.py').read_text()
(r/'payload.py').write_text(payload)
cmd=['ssh','-S','/Volumes/data/CodexTemp/T2/unfair-01-polynomials/macbook-control.sock','-o','BatchMode=yes','-o','StrictHostKeyChecking=yes','macbook-math','env PYTHONDONTWRITEBYTECODE=1 /Users/zhengzhihao/.venvs/math/bin/python -u -']
with (r/'stdout.jsonl').open('w') as out,(r/'stderr.txt').open('w') as err:
 p=subprocess.Popen(cmd,stdin=subprocess.PIPE,stdout=out,stderr=err,text=True)
 p.stdin.write(payload);p.stdin.close();rc=p.wait()
(r/'exit.json').write_text(json.dumps(dict(returncode=rc,payload_sha256=hashlib.sha256(payload.encode()).hexdigest(),input_sha256=hashlib.sha256((r/'inputs.json').read_bytes()).hexdigest(),command=cmd,finished_unix=time.time()),indent=2)+'\n')
raise SystemExit(rc)
