#!/usr/bin/env python3
"""Send and verify one small certificate in a single authenticated SSH session."""
import argparse
import base64
import hashlib
import json
from pathlib import Path
import subprocess
import time

p = argparse.ArgumentParser()
p.add_argument('certificate', type=Path)
p.add_argument('--verifier', type=Path, required=True)
p.add_argument('--out', type=Path, required=True)
a = p.parse_args()
payload = {}
for name, path in [('verify_matching_fast.py', a.verifier), ('certificate.json', a.certificate)]:
    raw = path.read_bytes()
    payload[name] = {'base64': base64.b64encode(raw).decode(), 'sha256': hashlib.sha256(raw).hexdigest()}
bootstrap = '''import base64,hashlib,json,subprocess,sys
from pathlib import Path
payload = json.loads(PAYLOAD)
root = Path('/Users/zhengzhihao/CodexRemoteMath/a28-selector/r005')
root.mkdir(parents=True,exist_ok=True)
for name, item in payload.items():
    raw = base64.b64decode(item['base64'])
    assert hashlib.sha256(raw).hexdigest()==item['sha256']
    (root/name).write_bytes(raw)
subprocess.run([sys.executable,str(root/'verify_matching_fast.py'),str(root/'certificate.json'),'--out',str(root/'remote-validation.json')],check=True,timeout=180)
'''.replace('PAYLOAD', repr(json.dumps(payload)))
begun = time.monotonic()
proc = subprocess.run(['ssh', '-o', 'BatchMode=yes', '-o', 'StrictHostKeyChecking=yes',
                       '-o', 'ConnectTimeout=20', 'macbook-math',
                       'nice -n 10 /usr/bin/python3 -'], input=bootstrap, text=True, capture_output=True)
result = {'ssh_exit_code': proc.returncode, 'elapsed_seconds': time.monotonic() - begun,
          'stdout': proc.stdout, 'stderr': proc.stderr,
          'input_sha256': {name: item['sha256'] for name, item in payload.items()}}
if proc.returncode == 0:
    result['verification'] = json.loads(proc.stdout.strip().splitlines()[-1])
    assert result['verification']['status'] == 'passed'
    assert result['verification']['certificate_sha256'] == payload['certificate.json']['sha256']
    assert result['verification']['verifier_sha256'] == payload['verify_matching_fast.py']['sha256']
a.out.write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result))
raise SystemExit(proc.returncode)
