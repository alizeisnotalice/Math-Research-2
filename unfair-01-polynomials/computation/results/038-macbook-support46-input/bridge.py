import base64,hashlib,json,os,subprocess,time,sys
from pathlib import Path
assert os.path.ismount('/Volumes/data')
r=Path(__file__).resolve().parent; smoke='--smoke' in sys.argv
out=r/('smoke' if smoke else 'run');out.mkdir(exist_ok=False)
s=(r/'frozen_support.py').read_text(); driver=(r/'remote.py').read_text()
payload='SOURCE_SHA='+repr(hashlib.sha256(s.encode()).hexdigest())+'\n'+s.split('\ndef main():')[0]+'\n'+driver+f'\nrun(smoke={smoke!r},workers={2 if smoke else 15})\n'
(out/'payload.py').write_text(payload)
cmd=['ssh','-o','BatchMode=yes','-o','StrictHostKeyChecking=yes','-o','ConnectTimeout=10','macbook-math','nice -n 10 env PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /Users/zhengzhihao/.venvs/math/bin/python -u -']
begin=time.time();count=0
with (out/'stderr.txt').open('w') as err, (out/'events.jsonl').open('w') as events:
 p=subprocess.Popen(cmd,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=err,text=True)
 p.stdin.write(payload);p.stdin.close()
 for line in p.stdout:
  row=json.loads(line)
  if row['event']=='chunk':
   raw=base64.b64decode(row.pop('gzip_base64'));assert hashlib.sha256(raw).hexdigest()==row['gzip_sha256']
   d=out/f"m{row['m']:02d}n{row['n']:02d}";d.mkdir(exist_ok=True)
   target=d/f"i{row['index']:06d}.jsonl.gz";assert not target.exists();target.write_bytes(raw)
   count+=1
  events.write(json.dumps(row,separators=(',',':'))+'\n');events.flush()
  if row['event']!='chunk' or count%250==0:print(json.dumps({k:v for k,v in row.items() if k!='unknown'}),flush=True)
 rc=p.wait()
(out/'exit.json').write_text(json.dumps(dict(returncode=rc,ssh_pid=p.pid,payload_sha256=hashlib.sha256(payload.encode()).hexdigest(),command=cmd,started_unix=begin,ended_unix=time.time(),chunks=count),indent=2)+'\n')
raise SystemExit(rc)
