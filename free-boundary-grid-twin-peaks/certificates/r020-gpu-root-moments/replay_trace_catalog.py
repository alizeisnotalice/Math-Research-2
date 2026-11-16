#!/usr/bin/env python3
"""Replay the frozen 62-case catalog into a caller-chosen temporary directory."""
import argparse, hashlib, json, os, subprocess, sys
from pathlib import Path
p=argparse.ArgumentParser()
p.add_argument('--binary',type=Path,required=True)
p.add_argument('--output-dir',type=Path,required=True)
p.add_argument('--limit',type=int,default=62)
a=p.parse_args(); base=Path(__file__).resolve().parent
assert 1<=a.limit<=62
out=a.output_dir.resolve(); out.mkdir(parents=True,exist_ok=True)
assert out != base, 'Use a temporary directory, not the certificate directory'
records=json.loads((base/'r020-frozen-62-certificate.json').read_text())['records'][:a.limit]
receipt=[]
for entry in records:
    cmd=entry['command'].copy()
    cmd[0]=str(a.binary.resolve()); cmd[6]=str(out/Path(cmd[6]).name)
    subprocess.run(cmd,check=True,cwd=out)
    trace=Path(cmd[6]+'.trace.json')
    subprocess.run([sys.executable,str(base/'verify_root_trace.py'),str(trace)],check=True,cwd=out)
    receipt.append({'case':entry['index'],'command':cmd,'trace_sha256':hashlib.sha256(trace.read_bytes()).hexdigest(),'passed':True})
(out/'replay-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'passed':len(receipt),'output_dir':str(out)}))
