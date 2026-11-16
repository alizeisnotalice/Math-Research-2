import json,sys,math,hashlib
from pathlib import Path
from collections import defaultdict
from decimal import Decimal as D,getcontext
getcontext().prec=80
outs=[]
for name in sys.argv[1:]:
 p=Path(name);raw=[json.loads(s,parse_float=D)for s in p.read_text().splitlines()];meta=raw[0];assert raw[-1]['kind']=='run_end'
 n=meta['n'];N=n*n;R=meta['replicates'];by=defaultdict(list);dom=int(meta['domain'],16)
 for x in raw[1:-1]:
  assert x['kind']=='replicate';t=x['tau'];k=x['replicate'];assert int(x['seed_hex'],16)==(0xA022000000000000+dom+n*1000000+t*10000+k)%(1<<64)
  if x['dead']:assert x['Z']==0 and x['logZ'] is None and x['U'] is None;x['U']=D(0)
  else:assert x['Z']>0 and x['U'] is not None
  by[t].append(x)
 assert sorted(by)==list(range(2,N+1))
 for t,v in by.items():assert sorted(x['replicate']for x in v)==list(range(R))
 Z=sum(sum(D(x['Z']) for x in v)/R for v in by.values());U=sum(sum(D(x['U'])for x in v)/R for v in by.values());Q=U/Z;var=D(0)
 for v in by.values():
  a=[D(x['U'])-Q*D(x['Z'])for x in v];mu=sum(a)/R;var+=sum((x-mu)**2 for x in a)/(R-1)/R
 expected_A=137077928808481869759712347602605837601908319657172480
 expected_D=2313395312707855007943501445611841392991226454185634816
 assert n==7;exactZ=D(expected_A)/math.factorial(N);exactQ=D(expected_D)/(N*expected_A)
 outs.append({'file':p.name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'metadata':meta,'records':sum(map(len,by.values())),'Z2':str(Z),'Q':str(Q),'stratified_delta_SE_diagnostic':str(var.sqrt()/Z),'exact_Q':str(exactQ),'Q_error':str(Q-exactQ),'Z2_relative_error':str(Z/exactZ-1),'dead_replicates':sum(x['dead']for v in by.values()for x in v),'wall_seconds':str(raw[-1]['wall_seconds']),'scope':'Finite calibration against exact n7; not a certified probability confidence interval or a mixing claim.'})
print(json.dumps(outs,indent=2,default=str))
