#!/usr/bin/env python3
"""Run exactly the authorized N128 area task and preserve its actual receipt."""
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
    ap=argparse.ArgumentParser();ap.add_argument('--producer',type=Path,required=True)
    ap.add_argument('--out-dir',type=Path,required=True);args=ap.parse_args();assert os.path.ismount('/Volumes/data')
    args.out_dir.mkdir(parents=True,exist_ok=True);frozen=args.out_dir/'sources/two_row_erosion_family.py'
    frozen.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(args.producer,frozen)
    cmd=[sys.executable,'-B',str(frozen),'--N','128','--bits','64','--out-dir',str(args.out_dir)]
    started=time.monotonic();p=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    samples=[];stop_reason=None
    while p.poll() is None:
        ps=subprocess.run(['ps','-p',str(p.pid),'-o','%cpu=,rss='],capture_output=True,text=True)
        if ps.returncode==0 and ps.stdout.strip():
            fields=ps.stdout.split();rss=int(fields[-1]);samples.append({'elapsed':time.monotonic()-started,'cpu_percent':float(fields[0]),'rss_kib':rss})
            if rss>512*1024:stop_reason='RSS_ABOVE_512MIB'
        if time.monotonic()-started>120:stop_reason='WALL_ABOVE_120SECONDS'
        if stop_reason:p.kill();break
        time.sleep(.1)
    stdout,stderr=p.communicate();gone=subprocess.run(['ps','-p',str(p.pid),'-o','pid='],capture_output=True,text=True)
    receipt={'command':cmd,'pid':p.pid,'exit_status':p.returncode,'wall_seconds':time.monotonic()-started,
             'owned_pid_absent':gone.returncode!=0 and not gone.stdout.strip(),'samples':samples,
             'stop_reason':stop_reason,'stdout':stdout,'stderr':stderr,'source_sha256':sha(frozen),
             'runner_sha256':sha(Path(__file__)),'python':sys.version,
             'scope':'ExactlyN128. NoN32rerun, no largerN, no MATLAB.'}
    result=args.out_dir/'result.json'
    if result.exists():
        d=json.loads(result.read_text());receipt['output_sha256']=sha(result);receipt['output_status']=d['status']
    (args.out_dir/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    assert p.returncode==0 and receipt['owned_pid_absent'] and receipt['output_status']=='PASS',receipt
    print(json.dumps({'status':'PASS','pid':p.pid,'seconds':receipt['wall_seconds'],
                      'sampled_peak_rss_kib':max((s['rss_kib'] for s in samples),default=None),
                      'result':str(result)}))


if __name__=='__main__':main()
