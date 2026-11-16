#!/usr/bin/env python3
"""Independent saved-mask and complete-interval audit; numerical evidence only."""
import json, hashlib, math, itertools, argparse, platform
from pathlib import Path
import numpy as np
ap=argparse.ArgumentParser();ap.add_argument('checkpoint',type=Path);ap.add_argument('roots',type=Path);ap.add_argument('source',type=Path);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--allow-partial',action='store_true');a=ap.parse_args()
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
c=json.loads(a.checkpoint.read_text());r=json.loads(a.roots.read_text())['candidates'][c['candidate_index']]
assert (c['status']=='complete' or a.allow_partial) and c['source_sha256']==sha(a.source) and c['roots_input_sha256']==sha(a.roots)
L=r['L'];blocks=r['blocks'];n=len(blocks);rows=c['E_progress'];assert sorted(e['E']for e in rows)==list(range(1,max(e['E']for e in rows)+1));assert a.allow_partial or len(rows)==L
C=np.array(r['C_coefficients_ascending']);assert set(C)<={0,1}
total=samples=survivors=compared_saved=missing_saved=0;maxres=0.;max_scaled=0.;seen=set();summaries=[]
for e in rows:
 E=e['E'];rad=([math.comb(L,E)]*n if E<L else [1]+[2]*(n-1));expected=math.prod(rad)
 assert e['status'] in ('complete','running','partial_unknown');assert (e['status']=='complete' or a.allow_partial) and e['radices_block0_first']==rad and e['total_codes']==expected
 hi=e['coverage_completed_half_open'][1];assert e['coverage_completed_half_open']==[0,hi] and 0<=hi<=expected and e['counts']['codes_seen']==hi;assert hi==expected if e['status']=='complete' else (a.allow_partial and E==max(z['E']for z in rows));assert len(e['top64'])==min(64,hi)
 assert 0<=e['counts']['closed_unit_interval_fractional_masks']<=e['counts']['fractional_masks']<=hi
 total+=hi
 options=sorted(itertools.combinations(range(L),E),key=lambda t:sum(1<<j for j in t)) if E<L else None
 local=[]
 for row in e['top64']:
  code=row['code'];assert 0<=code<hi and(E,code)not in seen;seen.add((E,code));v=code;digits=[];mask=0
  for j,b in enumerate(blocks):
   digit=v%rad[j];v//=rad[j];digits.append(digit)
   ids=([b['lift_groups'][z]for z in options[digit]] if E<L else ([b['original_group']]if digit==0 else b['lift_groups']))
   for z in ids:mask|=1<<z
  assert digits==row['mixed_radix_digits_block0_first'] and mask==int(row['A_root_group_bitmask_decimal'])
  A=np.array([1.]);B=np.array([1.])
  for j,g in enumerate(r['groups']):
   if mask&(1<<j):A=np.convolve(A,g['coefficients_ascending'])
   else:B=np.convolve(B,g['coefficients_ascending'])
  A/=A[0];B/=B[0];assert len(A)-1==row['degree_A'] and len(B)-1==row['degree_B']
  residual=float(abs(np.convolve(A,B)-C).max());maxres=max(maxres,residual)
  compared_saved+=int('cpu_float64_independent_root_reconstruction'in row);missing_saved+=int('cpu_float64_independent_root_reconstruction'not in row)
  saved=row.get('cpu_float64_independent_root_reconstruction',{'A_coefficients_ascending':A.tolist(),'B_coefficients_ascending':B.tolist()});err=max(float(abs(A-np.array(saved['A_coefficients_ascending'])).max()),float(abs(B-np.array(saved['B_coefficients_ascending'])).max()));scale=max(1.,float(abs(A).max()),float(abs(B).max()));max_scaled=max(max_scaled,err/scale);assert err/scale<1e-9
  values=np.r_[A,B];assert np.isfinite(values).all();mn=float(values.min());mx=float(values.max());fractional=bool(((values>1e-5)&(values<1-1e-5)).any());survivor=mn>=0 and mx<=1 and fractional;survivors+=survivor;samples+=1
  local.append({'code':code,'min':mn,'max':mx,'convolution_residual':residual,'numerical_fractional_survivor':survivor})
 summaries.append({'E':E,'codes':hi,'total_codes':expected,'gpu_counts':e['counts'],'samples':local})
full=sum(math.comb(L,E)**n for E in range(1,L))+2**(n-1);assert full==r['balanced_enumeration']['nonconstant_complement_representatives'];assert total<=full;assert a.allow_partial or total==full
out={'status':'partial_ledger_and_saved_samples_verified' if total<full else 'complete_ledger_and_saved_samples_verified','expected_full_codes':full,'candidate':r['name'],'codes':total,'saved_samples':samples,'saved_cpu_records_compared':compared_saved,'samples_without_saved_cpu_record':missing_saved,'saved_numerical_survivors':survivors,'max_numpy_convolution_residual':maxres,'max_scaled_saved_difference':max_scaled,'source_sha256':sha(a.source),'checkpoint_sha256':sha(a.checkpoint),'roots_sha256':sha(a.roots),'verifier_sha256':sha(Path(__file__)),'python':platform.python_version(),'numpy':np.__version__,'per_E':summaries,'claim_boundary':'Checks coverage bookkeeping and all saved masks independently. No uniform numerical error bound; no rigorous exclusion of unsaved masks.'}
a.out.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items()if k!='per_E'}))
