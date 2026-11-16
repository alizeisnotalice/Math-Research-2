#!/usr/bin/env python3
"""Small exact regression for general-n joint-null analytic formulas, not a proof by sampling."""
from pathlib import Path
from fractions import Fraction as F
import json,sys,platform,time,hashlib,os
p=Path(sys.argv[1]);cfg=json.loads(p.read_text());start=time.perf_counter();rows=[];prefix_count=radius_count=0
for n in range(cfg['minimum_n'],cfg['maximum_n']+1):
 e=F(1,n*n);h=F(n+1,n*n);K=n//2;m=K*K*e;old=(F(1,2*n)-e)**2;g=old+(1-2*e)*m/2
 gap=F(n-1,n**4)if n%2==0 else F(2*n-3,4*n**4)
 assert old==F((n-2)**2,4*n**4) and m/2-g==gap>0 and g<F(1,8)
 whole=3-3*e*e+4*g
 target=F(7,2)-F(4,n**3)+F(1,n**4)if n%2==0 else F(7,2)-F(1,n)+F(1,2*n*n)-F(2,n**3)
 assert whole==target<F(7,2)
 for j in range(2,n+1):
  margin=(j-1)*h+e-F(j,2*n);assert 2*n*n*margin==n*(j-2)+2*j>0;prefix_count+=1
 trunc=[]
 for mass in cfg['outside_masses']:
  b=F(mass);u2=min(F(1,4),(m+b)/4);v2=b/4
  assert 0<=v2<=u2<=F(1,4) and 0<=u2-v2<=m/4
  assert u2-v2==min((1-b)/4,m/4)
  trunc.append({'outside_mass':str(b),'u_squared':str(u2),'v_squared':str(v2),'increment':str(u2-v2)})
  radius_count+=1
 rows.append({'n':n,'K':K,'mP':str(m),'pure_internal_corner_area':str(old),'joint_null_corner_upper':str(g),'strict_half_mass_gap':str(gap),'complete_deletion_upper':str(whole),'truncated_radii':trunc})
assert next(r for r in rows if r['n']==3)['joint_null_corner_upper']=='5/108'
assert next(r for r in rows if r['n']==3)['complete_deletion_upper']=='85/27'
out={'schema':'r076-general-grid-joint-null-exact-scalar-regression-v1','passed':True,'n_count':len(rows),'actual_pure_prefix_inequality_checks':prefix_count,'actual_truncated_radius_checks':radius_count,'rows':rows,'scope':'Exact finite regression of parity/payment/pureprefix/truncatedradius identities only. General geometric theorem is analytic. No source-position search, jointprefix activity menu evaluation, actual input cap check, or halo area integration.','runtime':{'python':platform.python_version(),'platform':platform.platform(),'pid':os.getpid(),'seconds':time.perf_counter()-start,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'input_sha256':hashlib.sha256(p.read_bytes()).hexdigest()}}
Path(sys.argv[2]).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items()if k!='rows'}))
