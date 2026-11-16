from fractions import Fraction as F
from pathlib import Path
import json,sys,time,os,platform,hashlib
p=Path(sys.argv[1]);cfg=json.loads(p.read_text());start=time.time();out=[]
for case in cfg['cases']:
 n=case['n'];c=F(case['c']);e=F(1,n*n);h=F(n+1,n*n);q=n//2;side=F(1,2)-e;N=cfg['mesh'];inside=0;events=0
 for i in range(N):
  a=e+side*F(2*i+1,2*N)
  for j in range(N):
   b=e+side*F(2*j+1,2*N);accepted=False
   for d in sorted(set([a+k*h for k in range(q)]+[b+k*h for k in range(q)])):
    if d>=F(1,2):continue
    k=sum(ih*h<=d-a for ih in range(q));l=sum(ih*h<=d-b for ih in range(q))
    if not k*l:continue
    Lx=max(F(0),d-a+e);Ly=max(F(0),d-b+e)
    excess=F(k*l,n*n)-4*(1-c)*d*d-c*Lx*Ly;events+=1
    if excess>0:accepted=True;break
   inside+=accepted
 estimate=side*side*F(inside,N*N)
 out.append({'n':n,'candidate_density':str(c),'mesh':N,'midpoint_count':N*N,'accepted_midpoints':inside,'executed_event_tests':events,'midpoint_area_estimate':str(estimate),'area_estimate_decimal':float(estimate),'aux_corner_mass_over2':str(F(q*q,2*n*n)),'main_single_corner_threshold':str(F(1,8)+F(3,4*n**4))})
result={'schema':'r074-single-P-density-midpoint-probe-v1','cases':out,'scope':'Exact rational membership of sampled centers only; midpoint sums are NOT certified area bounds. Candidate densities are single-window upper limits, not yet all-square admissible. No finite atomic inference.','runtime':{'python':platform.python_version(),'pid':os.getpid(),'elapsed_seconds':time.time()-start,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'input_sha256':hashlib.sha256(p.read_bytes()).hexdigest()}}
Path(sys.argv[2]).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
