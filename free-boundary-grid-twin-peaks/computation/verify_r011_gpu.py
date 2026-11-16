#!/usr/bin/env python3
"""Recompute complete SIS run moments; numerical consistency, not certification."""
from pathlib import Path
from collections import defaultdict
import json,sys,math
B=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent/'certificates/r011-gpu-sis'
fields=['log_sumW','log_sumW2','log_sumWf','log_sumW2f','log_sumW2f2']
layerfields=['log_EW','log_EW2','log_EWf','log_EW2f','log_EW2f2']
seen=set();runs={};pools=defaultdict(list)
def calc(layers,n):
 by=defaultdict(lambda:[0.0]*6)
 for tau,M,mom in layers:
  d=by[tau];d[0]+=M
  for j,v in enumerate(mom):d[j+1]+=v
 a=sum(d[1]/d[0] for d in by.values());q=sum(d[3]/d[0] for d in by.values())/(a*n*n)
 var=0.;sq=0.
 for M,s,s2,sf,s2f,s2f2 in by.values():
  d2=s2f2/(n**4)-2*q*s2f/(n*n)+q*q*s2
  centered=d2-(sf/(n*n)-q*s)**2/M
  assert centered>=-1e-6*max(d2,1),(centered,d2)
  var+=max(0,centered)/(M*(M-1));sq+=s2/M**2
 return {'Q':q,'A2':a,'delta_se':math.sqrt(var)/a,'ESS':a*a/sq,'samples':int(sum(d[0] for d in by.values()))}
for n,count,batches in [(7,2,512),(8,8,768)]:
 for rep in range(count):
  name=f'n{n}-rep{rep:02d}-r11';rs=[json.loads(s) for s in (B/'runs'/name/'results.jsonl').read_text().splitlines()]
  meta=next(r for r in rs if r['kind']=='metadata');summary=next(r for r in rs if r['kind']=='summary')
  checks=[r for r in rs if r['kind']=='trajectory_check'];st=[r for r in rs if r['kind']=='stratum'];ba=[r for r in rs if r['kind']=='batch']
  assert len(checks)==len(st)==n*n-2 and all(r['ok'] for r in checks)
  assert len(ba)==batches and len({(r['tau'],r['index']) for r in ba})==batches
  assert {r['tau'] for r in st}==set(range(2,n*n))
  assert meta['thread_execution_width']==32 and meta['pilot_not_combined']
  layers=[]
  for s in st:
   rr=[r for r in ba if r['tau']==s['tau']];M=len(rr)*meta['batch'];assert M==s['samples']
   mom=[sum(math.exp(r[k]) if r[k] is not None else 0. for r in rr) for k in fields]
   for value,k in zip(mom,layerfields):
    if value:assert math.isclose(math.log(value/M),s[k],abs_tol=1e-10)
    else:assert s[k] is None
   assert sum(r['alive'] for r in rr)==s['alive']
   for r in rr:
    key=r['seed'];assert key not in seen;seen.add(key)
    assert r['bad_groups']==0 and 0<=r['alive']<=meta['batch']
    if r['alive']==0:assert all(r[k] is None for k in fields)
   layers.append((s['tau'],M,mom))
  ans=calc(layers,n);assert math.isclose(ans['Q'],summary['Q_euclidean_hat'],rel_tol=1e-11)
  assert math.isclose(math.log(ans['A2']),summary['log_A2_hat'],abs_tol=1e-10)
  pools[n].extend(layers);runs[name]=ans
pooled={str(n):calc(z,n) for n,z in pools.items()}
print(json.dumps({'status':'passed','scope':'numerical reaggregation and completeness; not an exact conditional probability','runs':runs,'pooled':pooled,'unique_batch_seeds':len(seen)},indent=2))
