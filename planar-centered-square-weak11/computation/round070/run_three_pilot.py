#!/usr/bin/env python3
"""Run exactly the three authorized R70 inputs with owned-process receipts."""
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


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--archive-root',type=Path,required=True)
    ap.add_argument('--spec',type=Path,required=True);ap.add_argument('--producer',type=Path,required=True)
    ap.add_argument('--out-dir',type=Path,required=True);a=ap.parse_args()
    assert os.path.ismount('/Volumes/data')
    sp=a.spec if a.spec.is_absolute() else a.archive_root/a.spec
    spec=json.loads(sp.read_text());assert [j['id'] for j in spec['jobs']]==['balanced8','uniform4-t099','historical16']
    a.out_dir.mkdir(parents=True,exist_ok=True);src=a.out_dir/'mass_sensitive_receiver.py'
    shutil.copyfile(a.producer,src);frozen_spec=a.out_dir/'three-pilot-spec.json';shutil.copyfile(sp,frozen_spec)
    def one(j):
        cmd=[sys.executable,'-B',str(src),'--archive-root',str(a.archive_root),'--spec',str(frozen_spec),
             '--job',j['id'],'--out-dir',str(a.out_dir)]
        start=time.perf_counter();p=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
        samples=[];reason=None
        while p.poll() is None:
            s=subprocess.run(['ps','-p',str(p.pid),'-o','rss='],capture_output=True,text=True)
            if s.returncode==0 and s.stdout.strip():
                rss=int(s.stdout.strip());samples.append({'elapsed':time.perf_counter()-start,'rss_kib':rss})
                if rss>1024*1024:reason='SAMPLED_RSS_OVER_1GIB'
            if time.perf_counter()-start>180:reason='WALL_OVER_180SECONDS'
            if reason:p.kill();break
            time.sleep(.1)
        out,err=p.communicate();gone=subprocess.run(['ps','-p',str(p.pid),'-o','pid='],capture_output=True,text=True)
        op=a.out_dir/(j['id']+'.json');cp=a.out_dir/(j['id']+'.cells.json.gz')
        rec={'id':j['id'],'command':cmd,'python':sys.version,'pid':p.pid,'exit_status':p.returncode,
             'owned_pid_absent':gone.returncode!=0 and not gone.stdout.strip(),'wall_seconds':time.perf_counter()-start,
             'stdout':out,'stderr':err,'stop_reason':reason,'samples':samples,'producer_sha256':sha(src),
             'spec_sha256':sha(frozen_spec),'input_sha256':j['sha256']}
        for k,q in [('output',op),('cells',cp)]:
            if q.exists():rec[k+'_sha256']=sha(q)
        (a.out_dir/(j['id']+'.receipt.json')).write_text(json.dumps(rec,indent=2)+'\n');return rec
    start=time.perf_counter()
    with ThreadPoolExecutor(max_workers=3) as pool:rows=list(pool.map(one,spec['jobs']))
    good=all(r['exit_status']==0 and r['owned_pid_absent'] and 'output_sha256' in r for r in rows)
    batch={'status':'PASS' if good else 'FAILED','runner_sha256':sha(Path(__file__)),'producer_sha256':sha(src),
           'spec_sha256':sha(frozen_spec),'python':sys.version,'wall_seconds':time.perf_counter()-start,
           'workers':3,'jobs':rows,'scope':'Exactly three directed R70 pilots; no repeated or random geometry jobs.'}
    (a.out_dir/'batch-receipt.json').write_text(json.dumps(batch,indent=2)+'\n')
    print(json.dumps({'status':batch['status'],'wall_seconds':batch['wall_seconds'],
                      'peak_sampled_rss_kib':max(s['rss_kib'] for r in rows for s in r['samples'])}))
    return 0 if good else 1


if __name__=='__main__':sys.exit(main())
