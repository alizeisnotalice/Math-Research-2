#!/usr/bin/env python3
"""Small exact replays plus authentication of the complete saved N9 scan."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time

parser=argparse.ArgumentParser()
parser.add_argument('--out',required=True,type=Path)
args=parser.parse_args()
root=Path(__file__).resolve().parent
out=args.out.resolve()
assert out!=root and root not in out.parents
out.mkdir(parents=True,exist_ok=False)
records=[]

def execute(name,cmd,cwd):
    start=time.monotonic()
    with (cwd/'stdout.txt').open('w') as stdout,(cwd/'stderr.txt').open('w') as stderr:
        child=subprocess.Popen(cmd,cwd=cwd,stdout=stdout,stderr=stderr)
        rc=child.wait()
    assert rc==0,(name,rc)
    records.append({'audit':name,'pid':child.pid,'wait_exit_code':rc,'wall_seconds':time.monotonic()-start})

specs=[
 ('dyadic-box','verify_dyadic_shell_samples.py','dyadic-shell-sample-check.json',False),
 ('dyadic-box','reinsertion_sum_small.py','reinsertion-sum-small.json',False),
 ('singleton-image','audit_singleton_image_n3.py','audit-singleton-image-n3.json',False),
 ('ring-sharpness','ring_gamma.py','ring_gamma.json',True),
]
for group,script,result,from_stdout in specs:
    work=out/script.removesuffix('.py');work.mkdir()
    shutil.copyfile(root/group/script,work/script)
    execute(script,[sys.executable,'-B',script],work)
    actual=json.loads((work/('stdout.txt' if from_stdout else result)).read_text())
    expected=json.loads((root/group/result).read_text())
    if group=='singleton-image':
        records[-1]['actual_python_field']=actual.pop('python')
        records[-1]['historical_python_field']=expected.pop('python')
    assert actual==expected,(group,script)
    records[-1]['exact_mathematical_json_match']=True
    records[-1]['source_sha256']=hashlib.sha256((work/script).read_bytes()).hexdigest()

nine=root/'nine-prefix'
work=out/'n9-evidence';work.mkdir()
execute('N9-full-saved-coverage',[sys.executable,'-B',str(nine/'scripts/audit_n9_certificate.py'),'--root',str(nine),'--json-out',str(work/'audit-result.json')],work)
assert json.loads((work/'audit-result.json').read_text())==json.loads((nine/'evidence/audit-result.json').read_text())
records[-1]['exact_mathematical_json_match']=True

work=out/'n9-min-executable';work.mkdir()
assert os.environ.get('CLANG_MODULE_CACHE_PATH'),'Set a shared compiler cache explicitly'
execute('N9-compiled-gates-and-16-hosts',['sh',str(nine/'scripts/min_replay.sh'),str(work/'run')],work)
receipt={'status':'passed','python':sys.version,'platform':platform.platform(),'scope':'Four small Python checks; full saved 326-checkpoint audit; compiled N7 gate, three N9 direct checks and 16 hosts only. No full N9 production rerun.','lifecycle':'Every child waited for and reaped. Ongoing project production jobs are outside this replay.','audits':records}
(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
