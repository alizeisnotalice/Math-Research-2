#!/usr/bin/env python3
"""Replay the small audits and saved n7 CRT in an external output directory."""
import argparse
import hashlib
import json
import pathlib
import platform
import shutil
import subprocess
import sys
import time

p=argparse.ArgumentParser()
p.add_argument('--out',type=pathlib.Path,required=True)
a=p.parse_args()
root=pathlib.Path(__file__).resolve().parent
out=a.out.resolve()
assert root not in out.parents and out!=root, 'Use an external output directory'
out.mkdir(parents=True,exist_ok=False)
records=[]

def run(group,script,result,stdout_result=False,remove_run=False):
    work=out/group
    work.mkdir()
    shutil.copyfile(root/group/script,work/script)
    begin=time.monotonic()
    with (work/'stdout.txt').open('w') as stdout, (work/'stderr.txt').open('w') as stderr:
        child=subprocess.Popen([sys.executable,'-B',script],cwd=work,stdout=stdout,stderr=stderr)
        rc=child.wait()
    assert rc==0,(group,rc)
    path=work/'stdout.txt' if stdout_result else work/result
    actual=json.loads(path.read_text())
    expected=json.loads((root/group/result).read_text())
    if remove_run:
        assert actual['run']['source_sha256']==expected['run']['source_sha256']
        actual.pop('run'); expected.pop('run')
    assert actual==expected,group
    records.append({'audit':group,'pid':child.pid,'wait_exit_code':rc,'wall_seconds':time.monotonic()-begin,'comparison':'exact mathematical JSON; run metadata excluded only for independent-path','source_sha256':hashlib.sha256((work/script).read_bytes()).hexdigest()})

run('path-encoding','audit_path_encoder_n3.py','path_encoder_n3_audit.json')
run('independent-path','audit_path_prefix.py','path-prefix-independent-r033.json',remove_run=True)
run('positive-moments','check_mu_coupling.py','mu-coupling-check.json')
run('two-seed','n3_pair_counts.py','n3_pair_counts.json',stdout_result=True)
# Both seed audits belong to one package but have separate work directories.
work=out/'two-seed-n4';work.mkdir()
begin=time.monotonic()
with (work/'stdout.txt').open('w') as stdout,(work/'stderr.txt').open('w') as stderr:
    child=subprocess.Popen([sys.executable,'-B',str(root/'two-seed/pair_completion_dp.py')],stdout=stdout,stderr=stderr)
    rc=child.wait()
assert rc==0
assert json.loads((work/'stdout.txt').read_text())==json.loads((root/'two-seed/n4_pair_completion.json').read_text())
records.append({'audit':'two-seed-n4','pid':child.pid,'wait_exit_code':rc,'wall_seconds':time.monotonic()-begin,'comparison':'exact JSON'})
n7=root/'n7-global-h45'
work=out/'n7';work.mkdir()
begin=time.monotonic()
cmd=[sys.executable,'-B',str(n7/'source/audit_n7_global_h45.py'),'--raw',str(n7/'raw'),'--source',str(n7/'source/short_gap_h45_mod.cpp'),'--pratt',str(n7/'references/pratt-certificates.json'),'--r31',str(n7/'references/r31-orbit-h1-h3-audit.json'),'--output',str(work/'audit.json')]
with (work/'stdout.txt').open('w') as stdout,(work/'stderr.txt').open('w') as stderr:
    child=subprocess.Popen(cmd,stdout=stdout,stderr=stderr);rc=child.wait()
assert rc==0
assert json.loads((work/'audit.json').read_text())==json.loads((n7/'reports/n7-global-h45-crt-audit.json').read_text())
records.append({'audit':'n7-saved-modular-CRT','pid':child.pid,'wait_exit_code':rc,'wall_seconds':time.monotonic()-begin,'comparison':'exact JSON; no production enumeration'})
receipt={'status':'passed','python':sys.version,'platform':platform.platform(),'method':'Subprocess.wait reaped every child; deterministic integer/Fraction audits','audits':records}
(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
