from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import subprocess

def sha(p):
    return sha256(Path(p).read_bytes()).hexdigest()

args=['/usr/bin/time','-l','./essential','5','8','6','8']
state=dict(command=args,source_sha256=sha('essential.cpp'),binary_sha256=sha('essential'),
           started_utc=datetime.now(timezone.utc).isoformat())
with open('n5-mac-b6-8.tsv','x') as out,open('n5-mac-b6-8.resource','x') as err:
    with subprocess.Popen(args,stdout=out,stderr=err) as p:
        state['pid']=p.pid
        Path('started.json').write_text(json.dumps(state,indent=2)+'\n')
        print('accepted',p.pid,flush=True)
        code=p.wait()
state.update(returncode=code,waited_reaped=True,ended_utc=datetime.now(timezone.utc).isoformat(),
             table_sha256=sha('n5-mac-b6-8.tsv'),resource_sha256=sha('n5-mac-b6-8.resource'))
Path('finished.json').write_text(json.dumps(state,indent=2)+'\n')
raise SystemExit(code)
