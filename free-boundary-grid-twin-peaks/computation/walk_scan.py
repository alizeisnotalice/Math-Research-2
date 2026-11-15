import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import subprocess
import time
from search import describe
p=argparse.ArgumentParser();p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--steps',type=int,default=5000000);p.add_argument('--workers',type=int,default=3);p.add_argument('--reflect-every',type=int,default=0);p.add_argument('--top-block-every',type=int,default=0);p.add_argument('--sizes',type=int,nargs='+',default=[4,5,8,12,16,24,32]);a=p.parse_args()
start=time.perf_counter()
jobs=[(n,k,kind) for n in a.sizes for k,kind in enumerate(('near','far'))]
def run(job):
 n,k,kind=job;seed=20261115000+10*n+k
 r=json.loads(subprocess.check_output([str(a.binary),str(n),str(a.steps),str(seed),kind,str(n*n),str(a.reflect_every),str(a.top_block_every)],text=True))
 cert=describe(r['best_labels'],n);assert cert['peak_count']==2 and cert['distance']==r['best_distance']
 r['best_verified']=True
 return r
out={'purpose':'conditional stationary walk for exploration; no asserted irreducibility, mixing, or rigorous confidence intervals',
 'binary_sha256':hashlib.sha256(a.binary.read_bytes()).hexdigest(),'workers':a.workers,'results':[]}
with ThreadPoolExecutor(max_workers=a.workers) as pool:
 for r in pool.map(run,jobs):
  out['results'].append(r);a.output.write_text(json.dumps(out,indent=2)+'\n')
  print(json.dumps({k:r[k] for k in ('n','start','accepted','samples','Q_empirical_unverified_mixing','best_distance')}),flush=True)
out['wall_seconds']=time.perf_counter()-start;a.output.write_text(json.dumps(out,indent=2)+'\n')
