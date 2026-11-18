"""Run a selected orbit range with exact output and completion hashes."""
from datetime import datetime,timezone
from hashlib import sha256
from pathlib import Path
import argparse
import json
import subprocess

D=Path(__file__).resolve().parent

def digest(p):
    return sha256(p.read_bytes()).hexdigest()

def now():
    return datetime.now(timezone.utc).isoformat()

def save(p,x):
    p.write_text(json.dumps(x,indent=2)+'\n')

p=argparse.ArgumentParser()
p.add_argument('name')
p.add_argument('lo',type=int)
p.add_argument('hi',type=int)
a=p.parse_args()
cmd=['/usr/bin/time','-l','./contact_geometry','5',str(a.hi-a.lo),str(a.lo),str(a.hi)]
raw=D/(a.name+'.tsv')
resource=D/(a.name+'.resource')
record=dict(command=cmd,started_utc=now(),source_sha256=digest(D/'contact_geometry.cpp'),
            binary_sha256=digest(D/'contact_geometry'),
            clang=subprocess.check_output(['clang++','--version'],text=True).splitlines()[0])
with raw.open('x') as output,resource.open('x') as errors:
    with subprocess.Popen(cmd,cwd=D,stdout=output,stderr=errors) as process:
        record['time_pid']=process.pid
        save(D/(a.name+'.started.json'),record)
        status=process.wait()
record.update(ended_utc=now(),waited_reaped=True,returncode=status,
              tsv_sha256=digest(raw),resource_sha256=digest(resource))
save(D/(a.name+'.finished.json'),record)
raise SystemExit(status)
