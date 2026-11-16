#!/usr/bin/env python3
"""Direct Decimal audit: complete tau strata, ratio of means, zero-inclusive delta."""
import csv,json,sys,hashlib,argparse
from decimal import Decimal as D,getcontext
from collections import defaultdict
from pathlib import Path
getcontext().prec=80

def audit(path,block_size=1):
 rows=list(csv.DictReader(Path(path).open())); n=int(rows[0]['n']);L=int(rows[0]['L']);by=defaultdict(list)
 R=max(int(x['rep']) for x in rows)+1
 domains={int(x['seed'])-int(x['tau'])*1000003-int(x['rep'])*1009-L*1000000007 for x in rows}
 assert len(domains)==1
 for x in rows:
  assert int(x['n'])==n and int(x['L'])==L
  by[int(x['tau'])].append(x)
 assert sorted(by)==list(range(2,n*n+1))
 for v in by.values():assert sorted(int(x['rep']) for x in v)==list(range(R))
 zm={t:sum(D(x['Z']) for x in v)/R for t,v in by.items()}
 um={t:sum(D(x['D'])/D(n*n) for x in v)/R for t,v in by.items()}
 Z=sum(zm.values());U=sum(um.values());Q=U/Z;var=D(0);ess={}
 for t,v in by.items():
  residuals=[D(x['D'])/D(n*n)-Q*D(x['Z']) for x in v];m=sum(residuals)/R
  var+=sum((u-m)**2 for u in residuals)/(R-1)/R
  zs=[D(x['Z']) for x in v];total=sum(zs)
  ess[t]={'ESS':str(total*total/sum(z*z for z in zs)) if total else '0','maxshare':str(max(zs)/total) if total else '0'}
 dead=sum(D(x['Z'])==0 for x in rows)
 assert all(D(x['D'])==0 for x in rows if D(x['Z'])==0)
 return {'file':Path(path).name,'sha256':hashlib.sha256(Path(path).read_bytes()).hexdigest(),'n':n,'L':L,'R':R,'records':len(rows),'seed_base':hex(next(iter(domains))),'Z2_count':str(Z),'U2_count_scaled':str(U),'Q_ratio_of_means':str(Q),'stratified_delta_SE_diagnostic':str(var.sqrt()/Z),'dead_replicates':dead,'zero_allowed_particle_events':sum(int(x['dead']) for x in rows),'block_size':block_size,'elementary_proposals':sum(int(x['proposals']) for x in rows),'accepted_blocks':sum(int(x['accept']) for x in rows),'block_acceptance':str(D(block_size)*sum(D(x['accept']) for x in rows)/sum(D(x['proposals']) for x in rows)) if L else None,'per_tau_weight_diagnostics':ess,'scope':'Finite numerical diagnostic; no certified confidence coverage or mixing claim. Main Q is ratio of stratum means, not the mean of replicate ratios.'}
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--block-size',type=int,default=1);ap.add_argument('files',nargs='+');a=ap.parse_args();assert a.block_size in (1,2)
 result={'status':'passed','arithmetic':'Decimal80 directly from CSV decimal strings','results':[audit(p,a.block_size) for p in a.files]};print(json.dumps(result,indent=2))
