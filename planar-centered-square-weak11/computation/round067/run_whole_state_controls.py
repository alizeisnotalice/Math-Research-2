#!/usr/bin/env python3
"""Exactly two authorized W/Wprime controls; no LP, theta or random expansion."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--state',type=Path,required=True)
    ap.add_argument('--module',type=Path,required=True);ap.add_argument('--erosion-source',type=Path,required=True)
    ap.add_argument('--padding-input',type=Path,required=True);ap.add_argument('--r64-input',type=Path,required=True)
    ap.add_argument('--out-dir',type=Path,required=True);args=ap.parse_args();assert os.path.ismount('/Volumes/data')
    assert sha(args.state)=='ab0e7dba2993869d69871384e9f95f0f095563b453dbde897561097d06251ccf'
    args.out_dir.mkdir(parents=True,exist_ok=True);frozen=args.out_dir/'sources/whole_erosion_state.py'
    frozen.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(args.state,frozen)
    records=[];started=time.monotonic()
    for name,ip in [('positive-grid-padding',args.padding_input),('r64-final',args.r64_input)]:
        op=args.out_dir/f'{name}.json';cmd=[sys.executable,'-B',str(frozen),'--module',str(args.module),
            '--erosion-source',str(args.erosion_source),'--input',str(ip),'--out',str(op)]
        t=time.monotonic();p=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
        try:stdout,stderr=p.communicate(timeout=60)
        except subprocess.TimeoutExpired:p.kill();p.communicate();raise
        ps=subprocess.run(['ps','-p',str(p.pid),'-o','pid='],capture_output=True,text=True)
        rec={'name':name,'command':cmd,'pid':p.pid,'exit_status':p.returncode,'wall_seconds':time.monotonic()-t,
             'stdout':stdout,'stderr':stderr,'owned_pid_absent':ps.returncode!=0 and not ps.stdout.strip(),
             'input_sha256':sha(ip),'state_source_sha256':sha(frozen)}
        records.append(rec)
        if op.exists():rec['output_sha256']=sha(op);rec['output_status']=json.loads(op.read_text())['status']
        assert p.returncode==0 and rec['owned_pid_absent'] and rec['output_status']=='PASS',rec
    receipt={'status':'PASS','runner_sha256':sha(Path(__file__)),'python':sys.version,'jobs':records,
             'wall_seconds':time.monotonic()-started,
             'scope':'Two directed prescribed inputs only; no additional state discount by Hremaining, no random inputs/LP/theta scan or MATLAB.'}
    (args.out_dir/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':'PASS','jobs':2,'seconds':receipt['wall_seconds']}))


if __name__=='__main__':main()
