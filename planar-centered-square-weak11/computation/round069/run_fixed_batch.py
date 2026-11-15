#!/usr/bin/env python3
"""Run the one authorized eight-job batch with per-process receipts."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def write(p,d):
    t=p.with_suffix(p.suffix+'.partial');t.write_text(json.dumps(d,indent=2)+'\n');t.replace(p)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--archive-root',type=Path,required=True)
    ap.add_argument('--spec',type=Path,required=True);ap.add_argument('--producer',type=Path,required=True)
    ap.add_argument('--out-dir',type=Path,required=True);ap.add_argument('--workers',type=int,default=8)
    args=ap.parse_args();assert os.path.ismount('/Volumes/data') and args.workers==8
    assert sha(args.producer)=='c7f64b45aa99f19e482e70479288c6e6e0bcaae0c780f35e983835630fba5ac7'
    sp=args.spec if args.spec.is_absolute() else args.archive_root/args.spec
    spec=json.loads(sp.read_text());assert len(spec['jobs'])==8
    args.out_dir.mkdir(parents=True,exist_ok=True);frozen=args.out_dir/'sources/central_peripheral_state.py'
    frozen.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(args.producer,frozen)
    fsp=args.out_dir/'spec-executed.json';shutil.copyfile(sp,fsp)
    def run(job):
        op=args.out_dir/f"{job['job_id']}.json"
        cmd=[sys.executable,'-B',str(frozen),'--archive-root',str(args.archive_root),
             '--spec',str(fsp),'--job',job['job_id'],'--out',str(op)]
        started=time.monotonic();p=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
        samples=[];stop=None
        while p.poll() is None:
            ps=subprocess.run(['ps','-p',str(p.pid),'-o','%cpu=,rss='],capture_output=True,text=True)
            if ps.returncode==0 and ps.stdout.strip():
                q=ps.stdout.split();rss=int(q[-1]);samples.append({'elapsed':time.monotonic()-started,'cpu_percent':float(q[0]),'rss_kib':rss})
                if rss>768*1024:stop='RSS_ABOVE_768MIB'
            if time.monotonic()-started>120:stop='WALL_ABOVE_120SECONDS'
            if stop:p.kill();break
            time.sleep(.05)
        stdout,stderr=p.communicate();ps=subprocess.run(['ps','-p',str(p.pid),'-o','pid='],capture_output=True,text=True)
        rec={'job_id':job['job_id'],'command':cmd,'pid':p.pid,'exit_status':p.returncode,
             'owned_pid_absent':ps.returncode!=0 and not ps.stdout.strip(),'wall_seconds':time.monotonic()-started,
             'source_sha256':sha(frozen),'spec_sha256':sha(fsp),'input_sha256':job['input_sha256'],
             'stdout':stdout,'stderr':stderr,'stop_reason':stop,'samples':samples}
        if op.exists():rec['output_sha256']=sha(op);rec['output_status']=json.loads(op.read_text())['status']
        write(args.out_dir/f"{job['job_id']}.receipt.json",rec);return rec
    started=time.monotonic()
    with ThreadPoolExecutor(max_workers=8) as pool:rows=list(pool.map(run,spec['jobs']))
    good=all(r['exit_status']==0 and r['owned_pid_absent'] and r.get('output_status')=='PASS' for r in rows)
    receipt={'status':'PASS' if good else 'FAILED','workers':8,'python':sys.version,'runner_sha256':sha(Path(__file__)),
             'source_sha256':sha(frozen),'spec_sha256':sha(fsp),'wall_seconds':time.monotonic()-started,'jobs':rows,
             'scope':'Exactly six historical maximum states plus two specified controls. No random samples, extra N, LP or MATLAB.'}
    write(args.out_dir/'batch-receipt.json',receipt)
    print(json.dumps({'status':receipt['status'],'jobs':8,'wall_seconds':receipt['wall_seconds'],
                      'sampled_peak_rss_kib':max((s['rss_kib'] for r in rows for s in r['samples']),default=None)}))
    return 0 if good else 1


if __name__=='__main__':sys.exit(main())
