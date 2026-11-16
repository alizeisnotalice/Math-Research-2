#!/usr/bin/env python3
"""Zero-inclusive ratio-of-stratum-means diagnostics from original-law CPU SMC."""
import csv,json,sys,hashlib
from pathlib import Path
from decimal import Decimal as D,getcontext
from collections import defaultdict
getcontext().prec=80
p=Path(sys.argv[1]);rows=list(csv.DictReader(p.open()));meta=json.loads(Path(str(p)+'.metadata.json').read_text())
n=int(meta['n']);N=n*n;R=int(meta['R']);L=int(meta['L']);by=defaultdict(list)
assert len(rows)==(N-1)*R
for x in rows:
 t=int(x['tau']);rep=int(x['rep']);assert int(x['n'])==n and int(x['L'])==L
 assert int(x['seed'])==int(meta['seed_base'])+t*1000003+rep*1009+L*1000000007
 by[t].append(x)
 z,d,f,g=[D(x[k])for k in ['Z','D','F','G']]
 assert z>=0 and d>=0 and f>=0
 assert abs(f-2*g-4*d)<=D('1e-12')*max(D(1),abs(f),abs(g),abs(d))
 assert abs(sum(D(x['COL'+str(j)])for j in range(n))-2*z)<=D('1e-12')*max(D(1),z)
 if not z:assert d==f==g==0
assert sorted(by)==list(range(2,N+1))
for t,v in by.items():assert sorted(int(x['rep'])for x in v)==list(range(R))
Z=sum(sum(D(x['Z']) for x in v)/R for v in by.values())
def ratio(key,scale):
 U=sum(sum(D(x[key])/scale for x in v)/R for v in by.values());q=U/Z;var=D(0)
 for v in by.values():
  es=[D(x[key])/scale-q*D(x['Z']) for x in v];m=sum(es)/R;var+=sum((e-m)**2 for e in es)/(R-1)/R
 return {'ratio':str(q),'stratified_delta_SE':str(var.sqrt()/Z)}
Q=ratio('D',D(N));V=ratio('F',D(8*N));C=ratio('G',D(4*N));assert abs(D(Q['ratio'])-2*(D(V['ratio'])-D(C['ratio'])))<D('1e-13')
col=[ratio('COL'+str(j),D(2))for j in range(n)];weights={}
for t,v in by.items():
 zz=[D(x['Z'])for x in v];s=sum(zz)
 weights[t]={'meanZ2share':str(s/R/Z),'ESS':str(s*s/sum(z*z for z in zz)) if s else '0','maxshare':str(max(zz)/s) if s else '0'}
out={'metadata':meta,'records':len(rows),'Z2_count':str(Z),'Q':Q,'V_over_n2':V,'C_over_n2':C,'column_marginal':col,'per_tau_weight_diagnostics':weights,'dead_replicates':sum(D(x['Z'])==0 for x in rows),'elementary_proposals':sum(int(x['proposals'])for x in rows),'csv_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'scope':'Finite double-precision particle estimate, aggregated by Decimal80. V,C use geometric center, whose population mean is exact by D4; no recentering by empirical mean. SE diagnostics are not certified confidence intervals or a mixing claim.'}
print(json.dumps(out,indent=2))
