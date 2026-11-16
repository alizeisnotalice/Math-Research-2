#!/usr/bin/env python3
"""Independent ledger/aggregate consistency checks, not asymptotic proof."""
import json, math
from pathlib import Path
from fractions import Fraction
from collections import Counter
B=Path(__file__).resolve().parent
C=B/'certificates'
def read(p): return json.loads(p.read_text())
def rows(p): return [json.loads(s) for s in p.read_text().splitlines()]
def close(a,b): assert math.isclose(a,b,rel_tol=1e-11,abs_tol=1e-12),(a,b)
out={'status':'passed','scope':'exact finite ledgers and numerical aggregate consistency only'}
A2={3:143808,4:779299174080,5:7063261481362958280000}
expected={3:Fraction(315,6848),4:Fraction(6063393111373,199500588564480),5:Fraction(403540943936240239007,41929158011105329152000)}
cut=C/'r010-natural-cuts';out['cuts']={}
for n in (3,4,5):
 d=read(cut/f'allsafe-n{n}.json');terms=d['by_size_cross']
 assert sum(int(x['raw_count']) for x in terms)==int(d['raw_sum'])
 w=sum((Fraction(int(x['raw_count']),2**x['cross_edges']) for x in terms),Fraction())
 assert w==Fraction(int(d['weighted_numerator']),int(d['common_denominator']))
 assert w/A2[n]==expected[n]
 out['cuts'][str(n)]=str(w/A2[n])
d=read(cut/'allsafe-n5-records.json');counts=Counter()
for r in d['records']:
 assert r['A_mask'].bit_count()==r['A_size']
 assert int(r['raw_count'])==math.comb(25,r['A_size'])*int(r['a_left'])*int(r['a_right'])
 counts[(r['A_size'],r['cross_edges'])]+=int(r['raw_count'])
assert counts=={(x['A_size'],x['cross_edges']):int(x['raw_count']) for x in read(cut/'allsafe-n5.json')['by_size_cross']}
assert len(d['records'])==10790
cpu=C/'r010-cpu-original-law';out['cpu']={}
for run in ('n7','n8-r001','n8-r002'):
 rs=rows(cpu/f'{run}-summary.jsonl');total=sum(float(r['A2_mean']) for r in rs)
 sq=sum(float(r['A2_mean'])**2/float(r['A2_ess']) for r in rs if float(r['A2_ess']))
 q=sum(float(r['A2_mean'])*float(r['Q_e2']) for r in rs)/total
 diag=read(cpu/f'{run}-global-diagnostics.json');agg=read(cpu/f'{run}-aggregate.json')
 close(total,diag['total_A2_estimate']);close(total**2/sq,diag['global_weight_ESS']);close(q,agg['Q_e2_ratio'])
 out['cpu'][run]={'Q':q,'ESS':total**2/sq,'samples':sum(r['samples'] for r in rs)}
saved=read(cpu/'saved-k2-labels.json');n=saved['n'];N=n*n
for r in saved['samples']:
 order=r['descending_order'];assert sorted(order)==list(range(N));lab=[0]*N
 for i,v in enumerate(order):lab[v]=N-i
 peaks=[]
 for v in range(N):
  x,y=v%n,v//n;adj=[Y*n+X for X,Y in ((x-1,y),(x+1,y),(x,y-1),(x,y+1)) if 0<=X<n and 0<=Y<n]
  if all(lab[v]>lab[u] for u in adj):peaks.append(v)
 assert len(peaks)==2 and set(peaks)==set(r['reported_peaks'])
 assert abs(peaks[0]%n-peaks[1]%n)+abs(peaks[0]//n-peaks[1]//n)==r['reported_rho']
out['independent_saved_label_checks']=len(saved['samples'])
gpu=C/'r010-gpu-original-law';out['gpu']={}
for n,nb in ((7,512),(8,768)):
 rs=rows(gpu/f'raw/adaptive-n{n}-r001/results.jsonl');st=[r for r in rs if r['kind']=='stratum'];ba=[r for r in rs if r['kind']=='batch'];su=next(r for r in rs if r['kind']=='summary')
 assert len(ba)==nb and {r['tau'] for r in st}==set(range(2,n*n))
 assert all(r['ok'] for r in rs if r['kind']=='trajectory_check')
 assert len({(r['tau'],r['index']) for r in ba})==nb
 assert all(r['alive']==0 for r in st if r['log_A2_tau_hat'] is None)
 a=max(r['log_A2_tau_hat'] for r in st if r['log_A2_tau_hat'] is not None);ws=[math.exp(r['log_A2_tau_hat']-a) if r['log_A2_tau_hat'] is not None else 0 for r in st]
 logz=a+math.log(sum(ws));q=sum(w*(r['Q_tau_euclidean'] or 0) for w,r in zip(ws,st))/sum(ws)
 close(logz,su['log_A2_hat']);close(q,su['Q_euclidean_hat'])
 assert sum(r['samples'] for r in st)==nb*65536
 out['gpu'][str(n)]={'Q':q,'samples':nb*65536}
rs=rows(C/'r010-gpu-kernel/production.jsonl');out['kernel']={}
for n in (7,8,9):
 rr=[r for r in rs if r.get('kind')=='replicate' and r['n']==n]
 assert len(rr)==512 and {r['replicate'] for r in rr}==set(range(512))
 a=max(r['logZ'] for r in rr);w=[math.exp(r['logZ']-a) for r in rr]
 out['kernel'][str(n)]={'scaled':sum(x*r['scaled'] for x,r in zip(w,rr))/sum(w),'Z_ESS':sum(w)**2/sum(x*x for x in w)}
assert all(r['ok'] for r in rs if r['kind']=='cpu_gpu_small_audit')
print(json.dumps(out,indent=2))
