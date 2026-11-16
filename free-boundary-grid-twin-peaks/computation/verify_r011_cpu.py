#!/usr/bin/env python3
"""50-digit reaggregation of floating CPU inputs; not exact probability data."""
import json,math
from pathlib import Path
from decimal import Decimal as D, getcontext
from collections import defaultdict
getcontext().prec=50
B=Path(__file__).resolve().parent/'certificates/r011-cpu-variance'
expected=['n3-domain0']+[f'n7-domain{i}' for i in range(4)]+[f'n8-domain{i}' for i in range(8)]
summary=json.loads((B/'all-variance-diagnostics.json').read_text());result={};pool=defaultdict(list);seen=set()
def aggregate(rows):
 layers=defaultdict(lambda:defaultdict(D))
 for r in rows:
  z=layers[r['tau']];z['N']+=r['samples']
  for k in ['sumW','sumW2','sumWf','sumW2f','sumW2f2']:z[k]+=D(r[k])
 A=sum(z['sumW']/z['N'] for z in layers.values());C=sum(z['sumWf']/z['N'] for z in layers.values());q=C/A
 V=D(0);sq=D(0)
 for z in layers.values():
  M=z['N'];d=z['sumWf']-q*z['sumW'];d2=z['sumW2f2']-2*q*z['sumW2f']+q*q*z['sumW2']
  variance_numerator=d2-d*d/M
  assert variance_numerator>=0,variance_numerator
  V+=variance_numerator/(M*(M-1));sq+=z['sumW2']/(M*M)
 return {'Q':q,'A2':A,'ESS':A*A/sq,'delta_se':V.sqrt()/A}
for name in expected:
 rs=[json.loads(x) for x in (B/name/'batches.jsonl').read_text().splitlines()];meta=rs.pop(0);n=meta['n'];pool[n].extend(rs)
 assert {r['tau'] for r in rs}==set(range(2,n*n+1))
 assert len(rs)==(n*n-1)*meta['batches']
 assert len({(r['tau'],r['batch']) for r in rs})==len(rs)
 for r in rs:
  assert r['samples']==meta['samples_per_tau_batch'] and r['raw_K2']+r['dead']==r['samples']
  assert r['seed'] not in seen;seen.add(r['seed'])
 a=aggregate(rs);s=summary['domains'][name]
 for k,field in [('Q','Q_e2'),('A2','A2_estimate'),('ESS','global_ess'),('delta_se','delta_method_se_Q')]:
  assert math.isclose(float(a[k]),s[field],rel_tol=1e-10,abs_tol=1e-12),(name,k,a[k],s[field])
 result[name]={k:str(v) for k,v in a.items()}
pooled={}
for n,rs in pool.items():
 a=aggregate(rs);assert math.isclose(float(a['Q']),summary['across_domains_by_n'][str(n)]['pooled_raw_moment_Q'],rel_tol=1e-12)
 pooled[str(n)]={k:str(v) for k,v in a.items()}
print(json.dumps({'status':'passed','scope':'higher-precision aggregation of previously rounded inputs; no confidence certificate','domains':result,'pooled':pooled,'unique_batch_seeds':len(seen)},indent=2))
