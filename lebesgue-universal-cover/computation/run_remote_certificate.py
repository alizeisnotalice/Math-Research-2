"""Stream the exact generator to the already configured SSH alias; no remote files.
Usage: python3 run_remote_certificate.py K rational_target output_json
"""
from pathlib import Path
import subprocess, sys, shlex, json, hashlib
source=Path(__file__).with_name('certify_five_body_upper.py').read_text()
k,target,output=sys.argv[1:]
runner="import sys,hashlib; s=sys.stdin.read(); exec(compile(s,'<streamed-certificate>','exec'),{'__name__':'__main__','STREAMED_SOURCE_SHA256':hashlib.sha256(s.encode()).hexdigest()})"
command='python3 -c '+shlex.quote(runner)+' '+shlex.quote(k)+' '+shlex.quote(target)
r=subprocess.run(['ssh','-o','BatchMode=yes','-o','StrictHostKeyChecking=yes','-o','ConnectTimeout=5','macbook-math',command],input=source,text=True,capture_output=True,timeout=600)
if r.returncode:raise RuntimeError((r.returncode,r.stderr))
x=json.loads(r.stdout)
assert x['script_sha256']==hashlib.sha256(source.encode()).hexdigest()
assert x['K']==int(k) and x['verified']
Path(output).write_text(r.stdout)
print(json.dumps({'exit_code':r.returncode,'K':x['K'],'strict_upper':x['strict_upper'],'display':x['strict_upper_display_only'],'points':x['sample_count'],'source_sha256':x['script_sha256']},indent=2))
