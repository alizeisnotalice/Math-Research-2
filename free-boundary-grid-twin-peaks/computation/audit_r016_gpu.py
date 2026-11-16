#!/usr/bin/env python3
"""Independent Decimal reaggregation of complete n10/n12 SMC records.
Empirical delta SE only; no mixing, PRNG independence or confidence certification.
"""
from pathlib import Path
from decimal import Decimal as D, getcontext
from collections import defaultdict
import json, hashlib, csv, sys, math
getcontext().prec=60
ROOT=Path(__file__).resolve().parents[1]
BASE=Path(sys.argv[1]) if len(sys.argv)>1 else ROOT/'certificates/r016-gpu'
ONE=Path(sys.argv[2]) if len(sys.argv)>2 else BASE/'onepeak'
seen=set()
def mv(x):
 n=len(x);assert n>1
 m=sum(x)/n
 return m,sum((v-m)**2 for v in x)/(n*(n-1))
def metal(p):return p.read_text().split('let metal = """',1)[1].split('"""',1)[0]
def load(p):
 rows=[json.loads(s,parse_float=D) for s in p.read_text().splitlines()]
 assert rows[0]['kind']=='metadata' and rows[-1]['kind']=='run_end'
 return rows[0],[r for r in rows if r['kind']=='replicate'],rows[-1]
def seed_check(m,x,one=False):
 expected=(int(m['seed_domain_hex'],16)+m['n']*1000000+x['replicate']+(0 if one else 0xA022000000000000+x['tau']*10000))%(1<<64)
 assert expected==int(x['seed_hex'],16) and expected not in seen
 seen.add(expected)
assert metal(BASE/'production_r013.swift')==metal(ONE/'source/onepeak_smc.swift')
assert hashlib.sha256(metal(BASE/'production_r013.swift').encode()).hexdigest()=='7474e22aee24049eb8692db789f32dff7b46e5c1c9d31f569622a7278c47da39'
strata=defaultdict(lambda:defaultdict(list));onevals=defaultdict(list);batches=[];telemetry=[]
for p in sorted((BASE/'production').glob('n*-p262144-r16-strat-domain*.jsonl')):
 m,reps,end=load(p); n=m['n'];R=m['replicates'];N=n*n
 assert m['particles']==262144 and m['tau_lo']==2 and m['tau_hi']==N
 groups=defaultdict(dict)
 for x in reps:
  seed_check(m,x);t=x['tau'];r=x['replicate'];assert r not in groups[t];groups[t][r]=x
  assert 0<=x['Z']<=1 and 0<=x['U']<=2*x['Z']
  if x['dead']:assert x['Z']==x['U']==0 and x['logZ'] is None
  else:assert x['Z']>0 and x['logZ'] is not None
 assert set(groups)==set(range(2,N+1))
 for t,g in groups.items():
  assert set(g)==set(range(R));strata[n][t].extend(g.values())
 z=sum(mv([x['Z'] for x in g.values()])[0] for g in groups.values())
 u=sum(mv([x['U'] for x in g.values()])[0] for g in groups.values());q=u/z
 vq=sum(mv([x['U']-q*x['Z'] for x in g.values()])[1] for g in groups.values())
 batches.append({'file':p.name,'n':n,'P':m['particles'],'R':R,'Z2':str(z),'Q':str(q),'empirical_delta_SE_Q':str(vq.sqrt()/z),'dead':sum(x['dead'] for x in reps),'nondead_zero_Z':sum(not x['dead'] and x['Z']==0 for x in reps),'wall_seconds':str(end['wall_seconds']),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
 mon=p.with_name(p.stem+'-monitor.csv')
 if mon.exists():
  rows=list(csv.DictReader(mon.open()));v=[D(x['agx_device_util_pct']) for x in rows if x.get('agx_device_util_pct')]
  telemetry.append({'file':mon.name,'samples':len(v),'mean_nonexclusive_device_percent':str(sum(v)/len(v)) if v else None,'samples_at_100':sum(x==100 for x in v)})
for p in sorted((ONE/'raw').glob('n*-p262144-r256-*.jsonl')):
 m,reps,end=load(p);n=m['n'];assert m['tau']==n*n+1 and m['particles']==262144
 assert [x['replicate'] for x in reps]==list(range(m['replicates']))
 for x in reps:
  seed_check(m,x,True);assert not x['dead'] and x['Z']>0 and x['tau']==n*n+1
 onevals[n].extend(x['Z'] for x in reps)
assert set(strata)==set(onevals)=={10,12}
pooled=[]
for n,g in sorted(strata.items()):
 assert all(len(v)==32 for v in g.values()) and len(onevals[n])==512
 z2=sum(mv([x['Z'] for x in v])[0] for v in g.values());u2=sum(mv([x['U'] for x in v])[0] for v in g.values());q=u2/z2
 varz=sum(mv([x['Z'] for x in v])[1] for v in g.values())
 varq=sum(mv([x['U']-q*x['Z'] for x in v])[1] for v in g.values())
 z1,var1=mv(onevals[n]);ratio=z2/z1
 vr=varz/z1**2+z2**2*var1/z1**4
 pooled.append({'n':n,'P':262144,'twopeak_repeats_per_tau':32,'onepeak_repeats':512,'Z1':str(z1),'Z2':str(z2),'Q':str(q),'Q_empirical_delta_SE':str(varq.sqrt()/z2),'R':str(ratio),'R_empirical_delta_SE':str(vr.sqrt()),'R_over_N':str(ratio/D(n*n)),'logR_over_n':str(ratio.ln()/D(n)),'logR_over_nlogn':str(ratio.ln()/(D(n)*D(n).ln()))})
cal_m,cal_reps,cal_end=load(ONE/'raw/n7-calibration.jsonl')
assert cal_m['n']==7 and cal_m['tau']==50 and cal_m['particles']==8192
assert [x['replicate'] for x in cal_reps]==list(range(16))
for x in cal_reps:
 seed_check(cal_m,x,True);assert not x['dead'] and x['Z']>0 and x['tau']==50
cal_mean,cal_var=mv([x['Z'] for x in cal_reps])
cal_exact=D(134267527713152677457983522684235911142268448963584)/D(math.factorial(49))
calibration={'n':7,'P':8192,'R':16,'Z1':str(cal_mean),'exact_Z1':str(cal_exact),'relative_difference':str(cal_mean/cal_exact-1)}
out={'status':'passed','calibration':calibration,'decimal_precision':60,'unique_seed_records':len(seen),'canonical_Metal_hash':hashlib.sha256(metal(BASE/'production_r013.swift').encode()).hexdigest(),'batches':batches,'pooled':pooled,'telemetry':telemetry,'evidence':'Finite numerical evidence. SE descriptive, not confidence bounds. Seed uniqueness does not prove independence.'}
print(json.dumps(out,indent=2))
