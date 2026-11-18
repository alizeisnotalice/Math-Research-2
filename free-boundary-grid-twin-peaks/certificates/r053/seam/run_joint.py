"""Run the three bounded modular rectangle workers and preserve receipts."""
from argparse import ArgumentParser
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import subprocess

HERE=Path(__file__).resolve().parent

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def run(ix,folder,plan):
    out=folder/f'prime-{ix}.json';log=folder/f'prime-{ix}.resource'
    argv=['/usr/bin/time','-l',str(HERE/'joint_seam'),'6','7',str(ix),'--max-table-mib','40']
    with out.open('x') as f,log.open('x') as g:
        p=subprocess.Popen(argv,stdout=f,stderr=g,start_new_session=True)
        try:
            code=p.wait(timeout=40)
        except subprocess.TimeoutExpired:
            os.killpg(p.pid,signal.SIGTERM)
            try:p.wait(timeout=2)
            except subprocess.TimeoutExpired:
                os.killpg(p.pid,signal.SIGKILL)
                p.wait()
            raise
    assert code==0,log
    r=json.loads(out.read_text())
    assert r['mode']=='count' and r['prime_index']==ix and (r['width'],r['height'])==(6,7)
    assert r['processed_by_layer']==plan['states_by_layer'] and r['max_table_bytes']==plan['max_table_bytes']
    rss=re.search(r'^\s*(\d+)\s+maximum resident set size$',log.read_text(),re.M)
    assert rss
    return {'prime_index':ix,'waited':True,'exit_status':code,'pid':p.pid,'max_rss_bytes':int(rss[1]),
            'stdout_sha256':digest(out),'resource_sha256':digest(log),'command':argv}


def main():
    a=ArgumentParser();a.add_argument('output',type=Path);b=a.parse_args()
    b.output.mkdir(exist_ok=False)
    plan=json.loads((HERE/'n12-plan.json').read_text())
    assert plan['independent_sets']==5_598_861 and plan['max_table_bytes']==38_252_064
    with ThreadPoolExecutor(max_workers=3) as pool:
        result=list(pool.map(lambda i:run(i,b.output,plan),range(3)))
    receipt={'counter_source_sha256':digest(HERE/'joint_seam.cpp'),'binary_sha256':digest(HERE/'joint_seam'),
             'plan_sha256':digest(HERE/'n12-plan.json'),'jobs':result,'max_worker_rss_bytes':max(r['max_rss_bytes'] for r in result)}
    (b.output/'receipts.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print('Waited all modular rectangle workers, max worker RSS',receipt['max_worker_rss_bytes'])

if __name__=='__main__':main()
