#!/usr/bin/env python3
import hashlib,json,math,sys
from pathlib import Path
BASE='934275c0da831ea81c22ac78fcd4175517723afef9343553796f77e71d3bce63'
SOURCE='3a9f3c0a8c10409346c25f16f2d82463f928bfaacf7307b4a88a770ec8f69ec6'
TOP=Path(__file__).resolve().parent/'exact-top2-d4'
def ess(w): return sum(w)**2/sum(x*x for x in w)
def main():
 if len(sys.argv) not in (3,4): raise SystemExit('usage: aggregate_calibration.py raw.jsonl source.swift [exact-top2-dir]')
 raw,src=map(Path,sys.argv[1:3]); top=Path(sys.argv[3]) if len(sys.argv)==4 else TOP
 lines=[json.loads(x) for x in raw.read_text().splitlines()]; meta=lines[0]
 assert hashlib.sha256(src.read_bytes()).hexdigest()==SOURCE and meta['implementation']=='r026-tau2-firstmerge-rb-terminal'
 for k in ('n','particles','replicates','tau_lo','tau_hi','move_steps','domain','method'): assert k in meta
 assert meta['tau_lo']==meta['tau_hi']==2 and meta['method']=='strat' and meta['move_steps']==16
 assert lines[-1]['kind']=='run_end' and math.isfinite(lines[-1]['wall_seconds'])
 reps={x['replicate']:x for x in lines if x['kind']=='replicate'}; rb={x['replicate']:x for x in lines if x['kind']=='firstmerge_rb_calibration'}
 assert set(reps)==set(rb)==set(range(meta['replicates']))
 n=meta['n']; dom=int(meta['domain'],16)
 for r,x in reps.items():
  assert x['tau']==2 and not x['dead'] and math.isfinite(x['Z']) and x['Z']>0 and math.isfinite(x['logZ'])
  seed=(0xA022000000000000+dom+n*1000000+20000+r)&((1<<64)-1); assert int(x['seed_hex'],16)==seed
  y=rb[r]; h=y['mean_m_hist']; assert len(h)==2*n-1 and all(math.isfinite(v) and v>=0 for v in h)
  assert all(math.isfinite(y[k]) for k in ('mean_m','mean_m_d2','mean_m_F0','mean_m_G0')) and abs(sum(h)-y['mean_m'])<1e-11
  assert abs(y['mean_m_F0']-2*y['mean_m_G0']-4*y['mean_m_d2'])<1e-10
 reps=[reps[i] for i in range(meta['replicates'])]; rb=[rb[i] for i in range(meta['replicates'])]
 mx=max(x['logZ'] for x in reps); w=[math.exp(x['logZ']-mx) for x in reps]; s=sum(w)
 av=lambda k:sum(a*b[k] for a,b in zip(w,rb))/s
 assert n in (3,7), 'no silent fallback exact baseline'
 exact_path=top/f'n{n}.json'; ex=json.loads(exact_path.read_text()); assert ex['n']==n
 A=int(ex['A2']); B=int(ex['Btop2']); assert B>0 and A>0
 summary=json.loads((top/'summary-n3-n7.json').read_text())
 cross=next(x for x in summary['results'] if x['n']==n); assert int(cross['Btop2'])==B and int(cross['K2_A2'])==A
 ED=int(ex['squared_euclidean_total'])/A; EF=int(ex['firstsum'])/A; EG=int(ex['secondsum'])/A
 hist=[int(ex['distance_counts'].get(str(d),0))/A for d in range(2*n-1)]
 got=[sum(w[i]*rb[i]['mean_m_hist'][d] for i in range(len(w)))/s for d in range(2*n-1)]
 Bhat=math.exp(mx)*s/len(w); Ahat=Bhat*av('mean_m'); fact=math.factorial(n*n)
 out={'audit_passed':True,'audit_scope':'input/source/count consistency only, not numerical accuracy','actual_calibration_source_sha256':SOURCE,'frozen_trace_source_sha256':BASE,'raw_sha256':hashlib.sha256(raw.read_bytes()).hexdigest(),'exact_source_file':exact_path.name,'exact_source_sha256':hashlib.sha256(exact_path.read_bytes()).hexdigest(),'summary_cross_source_sha256':hashlib.sha256((top/'summary-n3-n7.json').read_bytes()).hexdigest(),'n':n,'P':meta['particles'],'R':len(reps),'logZ_weight_ESS':ess(w),'logZ_weight_maxshare':max(w)/s,'B_over_Nfactorial_estimate':Bhat,'mA2_over_Nfactorial_estimate':Ahat,'mD_over_Nfactorial_estimate':Bhat*av('mean_m_d2'),'mF_over_Nfactorial_estimate':Bhat*av('mean_m_F0'),'mG_over_Nfactorial_estimate':Bhat*av('mean_m_G0'),'B_over_Nfactorial_exact':B/fact,'A2_over_Nfactorial_exact':A/fact,'B_relative_error':Bhat/(B/fact)-1,'A2_relative_error':Ahat/(A/fact)-1,'A2_over_B_exact':A/B,'A2_over_B_mRB_estimate':av('mean_m'),'Q_mRB':av('mean_m_d2')/(n*n*av('mean_m')),'Q_L_legacy':sum(a*b['tau2_lift_scaled_mean_L_distance2'] for a,b in zip(w,reps))/s/(n*n*sum(a*b['tau2_lift_scaled_mean_L'] for a,b in zip(w,reps))/s),'exact':{'A2':str(A),'Btop2':str(B),'Q':ED/(n*n),'D':ED,'F':EF,'G':EG,'distance':hist},'mRB':{'D':av('mean_m_d2')/av('mean_m'),'F':av('mean_m_F0')/av('mean_m'),'G':av('mean_m_G0')/av('mean_m'),'distance':[x/av('mean_m') for x in got]},'note':'finite SMC diagnostic only; no SE pass criterion and no cone claim.'}
 out['differences']={k:out['mRB'][k]-out['exact'][k] for k in ('D','F','G')}
 print(json.dumps(out,indent=2))
if __name__=='__main__': main()
