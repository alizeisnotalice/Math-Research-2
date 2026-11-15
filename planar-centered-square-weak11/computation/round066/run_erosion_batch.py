#!/usr/bin/env python3
"""Bounded authorized erosion jobs with actual process receipts and RSS samples."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def write(p,d):
    p.parent.mkdir(parents=True,exist_ok=True);t=p.with_suffix(p.suffix+'.partial')
    t.write_text(json.dumps(d,indent=2)+'\n');t.replace(p)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--budget-root',type=Path,required=True)
    ap.add_argument('--module',type=Path,required=True);ap.add_argument('--erosion',type=Path,required=True)
    ap.add_argument('--out-dir',type=Path,required=True);ap.add_argument('--phase',choices=('pilot','strata'),required=True)
    ap.add_argument('--workers',type=int,default=12);ap.add_argument('--seconds-per-job',type=int,default=120)
    ap.add_argument('--rss-limit-mib',type=int,default=768);args=ap.parse_args()
    assert os.path.ismount('/Volumes/data') and 1<=args.workers<=12
    assert sha(args.erosion)=='4dc4fb494736c40334359f9950f052ee2e4a53b503883ffbc3f9bc0da69b1361'
    assert sha(args.module)=='4ef61b815b9c4a9a8972cdaf31716ae7917b59243a29fb8697c8042137788246'
    args.out_dir.mkdir(parents=True,exist_ok=True);frozen=args.out_dir/'sources/erosion_ring_diagnostic.py'
    frozen.parent.mkdir(parents=True,exist_ok=True)
    if frozen.exists():assert sha(frozen)==sha(args.erosion)
    else:shutil.copyfile(args.erosion,frozen)
    jobs=[('r64-final',args.budget_root/'r64-final-theta.json')] if args.phase=='pilot' else [(f'stratum-{i:02}',args.budget_root/f'stratum-{i:02}.json') for i in range(12)]
    def run(job):
        name,source=job;out=args.out_dir/f'{name}-erosion.json'
        cmd=[sys.executable,'-B',str(frozen),'--module',str(args.module),'--budget-result',str(source),'--out',str(out)]
        start=time.monotonic();proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
        receipt={'job':name,'command':cmd,'pid':proc.pid,'source_sha256':sha(frozen),
                 'geometry_sha256':sha(args.module),'budget_result_sha256':sha(source),'samples':[]}
        stop_reason=None
        while proc.poll() is None:
            ps=subprocess.run(['ps','-p',str(proc.pid),'-o','%cpu=,rss='],capture_output=True,text=True)
            if ps.returncode==0 and ps.stdout.strip():
                parts=ps.stdout.split();rss=int(parts[-1]);receipt['samples'].append({'elapsed_seconds':time.monotonic()-start,'cpu_percent':float(parts[0]),'rss_kib':rss})
                if rss>args.rss_limit_mib*1024:stop_reason='RSS_LIMIT_EXCEEDED'
            if time.monotonic()-start>args.seconds_per_job:stop_reason='WALL_LIMIT_EXCEEDED'
            if stop_reason:proc.kill();break
            time.sleep(.05)
        stdout,stderr=proc.communicate();absent=subprocess.run(['ps','-p',str(proc.pid),'-o','pid='],capture_output=True,text=True)
        receipt.update({'exit_status':proc.returncode,'wall_seconds':time.monotonic()-start,'stdout':stdout,'stderr':stderr,
                        'owned_pid_absent':absent.returncode!=0 and not absent.stdout.strip(),'stop_reason':stop_reason})
        if out.exists():
            result=json.loads(out.read_text());receipt.update({'output_sha256':sha(out),'output_status':result['status'],
                    'cases':len(result['cases']),'blocks':sum(len(c['blocks']) for c in result['cases'])})
        write(args.out_dir/f'{name}.receipt.json',receipt);return receipt
    started=time.monotonic()
    with ThreadPoolExecutor(max_workers=args.workers) as pool:receipts=list(pool.map(run,jobs))
    payload={'status':'PASS' if all(r['exit_status']==0 and r['owned_pid_absent'] and r.get('output_status')=='PASS' for r in receipts) else 'FAILED',
             'phase':args.phase,'workers':min(args.workers,len(jobs)),'wall_seconds':time.monotonic()-started,
             'runner_sha256':sha(Path(__file__)),'python':sys.version,'platform':platform.platform(),
             'limits':{'seconds_per_job':args.seconds_per_job,'rss_limit_mib':args.rss_limit_mib},'jobs':receipts}
    write(args.out_dir/f'batch-{args.phase}.receipt.json',payload)
    print(json.dumps({'status':payload['status'],'phase':args.phase,'workers':payload['workers'],'seconds':payload['wall_seconds'],'cases':sum(r.get('cases',0) for r in receipts),'blocks':sum(r.get('blocks',0) for r in receipts)}))
    if payload['status']!='PASS':return 1
    return 0


if __name__=='__main__':sys.exit(main())
