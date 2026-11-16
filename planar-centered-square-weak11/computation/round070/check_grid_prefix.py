#!/usr/bin/env python3
from pathlib import Path
from fractions import Fraction as F
from bisect import bisect_left,bisect_right
from concurrent.futures import ProcessPoolExecutor
import os,json,math,time,hashlib,sys,subprocess

def job(arg):
 n, nums, stripe, stripes=arg;k=math.isqrt(n);D=32*k*n*n;R=(k-1)*D//(2*k)
 coords=[j*(n+1)*D//(n*n) for j in range(n)];rows=[];all_events=0;t=time.perf_counter()
 def count(X,c,strict=False):
  return (bisect_left(coords,X+c)-bisect_right(coords,X-c)) if strict else (bisect_right(coords,X+c)-bisect_left(coords,X-c))
 for ix,xnum in enumerate(nums):
  if ix%stripes!=stripe:continue
  X=xnum*D//16
  for ynum in nums:
   Y=ynum*D//16;U=count(X,R,True)*count(Y,R,True)
   active=0;failed=0;first=None;mn=None;events=sorted({abs(z-X) for z in coords}|{abs(z-Y) for z in coords})
   feasible=0<U<n*n and U*k*k<(k-1)**2*n*n
   for C in events:
    if C<R:continue
    H=count(X,C)*count(Y,C)-U
    assert H>=0
    if H==0:continue
    all_events+=1
    union_num=(k-1)**2*n*n+H*k*k
    if union_num<=k*k*n*n:continue
    active+=1
    slack=(C+R)**2*k*k*n*n-D*D*union_num
    mn=slack if mn is None else min(mn,slack)
    if slack<=0:
     failed+=1
     if first is None:first={'distance':str(F(C,D)),'prefix_mass':str(F(H,n*n)),'squared_slack':str(F(slack,D*D*k*k*n*n))}
   rows.append({'x':str(F(xnum,16)),'y':str(F(ynum,16)),'U_count':U,'candidate_mass_exceeds_U':feasible,'active_prefix_count':active,'violated_prefix_count':failed,'all_prefix_pass':feasible and failed==0,'min_active_squared_slack':None if mn is None else str(F(mn,D*D*k*k*n*n)),'first_violation':first})
 return {'n':n,'stripe':stripe,'pid':os.getpid(),'events_checked':all_events,'rows':rows,'elapsed_seconds':time.perf_counter()-t}

def main():
 assert os.path.ismount('/Volumes/data') and __debug__
 p=Path(__file__).resolve().parent;inp=p/'r070-grid-prefix-input.json';cfg=json.loads(inp.read_text());start=time.perf_counter()
 jobs=[(n,cfg['coordinate_numerators'],j,cfg['stripes']) for n in cfg['grid_sizes'] for j in range(cfg['stripes'])]
 with ProcessPoolExecutor(max_workers=cfg['workers']) as ex:results=list(ex.map(job,jobs))
 pids=sorted({r['pid'] for r in results});gone={str(pid):subprocess.run(['ps','-p',str(pid),'-o','pid='],capture_output=True,text=True).returncode!=0 for pid in pids}
 assert all(gone.values())
 summary=[]
 for n in cfg['grid_sizes']:
  rr=[z for r in results if r['n']==n for z in r['rows']]
  assert len(rr)==len(cfg['coordinate_numerators'])**2
  summary.append({'n':n,'sources':n*n,'candidate_m':str(F((math.isqrt(n)-1)**2,math.isqrt(n)**2)),'point_count':len(rr),'all_prefix_pass_count':sum(z['all_prefix_pass'] for z in rr),'prefix_violation_point_count':sum(z['violated_prefix_count']>0 for z in rr),'events_checked':sum(r['events_checked'] for r in results if r['n']==n)})
 result={'schema':'r070-grid-all-prefix-point-check-v1','scope':cfg['scope'],'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'input_sha256':hashlib.sha256(inp.read_bytes()).hexdigest(),'python':sys.version,'pid':os.getpid(),'workers_requested':cfg['workers'],'worker_pids':pids,'all_owned_worker_pids_absent':gone,'wall_seconds':time.perf_counter()-start,'summary':summary,'chunks':results}
 out=p/'r070-grid-prefix-result.json';out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k] for k in ['summary','worker_pids','wall_seconds']}),flush=True)
if __name__=='__main__':main()
