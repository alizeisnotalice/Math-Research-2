#!/usr/bin/env python3
"""Throughput-oriented v2 of the root-group screen; v1 remains untouched.

This imports only schema/CPU-recheck helpers from v1.  The batch recurrence is
implemented here so that MPS allocations are not emptied per batch and top-k
metrics are copied to CPU as one packed tensor rather than per candidate.
"""
from __future__ import annotations
import argparse, hashlib, importlib.util, json, platform, sys
from pathlib import Path
from time import perf_counter
import numpy as np
import torch

V1 = Path(__file__).with_name('gpu_root_group_search.py')
spec = importlib.util.spec_from_file_location('root_group_v1', V1); v1 = importlib.util.module_from_spec(spec)
assert spec and spec.loader; spec.loader.exec_module(v1)

def atomic(path, value):
    tmp=path.with_suffix(path.suffix+'.tmp'); tmp.write_text(json.dumps(value,indent=2)+'\n'); tmp.replace(path)

def tmul(poly, factor):
    out=torch.zeros_like(poly); width=poly.shape[1]
    for shift in range(factor.numel()): out[:,shift:]+=factor[shift]*poly[:,:width-shift]
    return out

def load_v2(path: Path, candidate_index: int):
    """Prepared-container loader with v2's explicit 30-group ceiling."""
    raw=json.loads(path.read_text()); candidates=raw.get('candidates')
    if not isinstance(candidates,list) or not 0<=candidate_index<len(candidates): raise SystemExit('invalid candidate index')
    candidate=candidates[candidate_index]; C=np.asarray(candidate.get('coefficients_ascending'),dtype=np.float64); factors=candidate.get('groups_ascending_coefficients')
    if C.ndim!=1 or len(C)<2 or not isinstance(factors,list) or not 2<=len(factors)<=30: raise SystemExit('invalid candidate or group_count outside v2 limit 2..30')
    fs=[np.asarray(x,dtype=np.float64) for x in factors]
    if any(f.ndim!=1 or len(f) not in (2,3) or f[-1]!=1.0 for f in fs): raise SystemExit('groups must be monic real degree-1/2 ascending factors')
    product=np.array([1.0])
    for f in fs: product=np.convolve(product,f)
    err=np.max(np.abs(product-C))/max(1.0,np.max(np.abs(C)))
    if err>1e-8 or abs(C[0])<1e-14: raise SystemExit(f'prepared product mismatch ({err}) or zero C(0)')
    return raw,candidate,C,fs,err

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('input',type=Path); ap.add_argument('--candidate-index',type=int,required=True); ap.add_argument('--checkpoint',type=Path,required=True)
    ap.add_argument('--start',type=int,default=0); ap.add_argument('--end',type=int); ap.add_argument('--batch-size',type=int,default=131072)
    ap.add_argument('--keep',type=int,default=64); ap.add_argument('--device',choices=('cpu','mps'),default='cpu'); ap.add_argument('--tol',type=float,default=1e-6); ap.add_argument('--eta',type=float,default=1e-5); ap.add_argument('--condition-limit',type=float,default=1e8); ap.add_argument('--smoke',action='store_true')
    a=ap.parse_args()
    if a.batch_size<1 or a.keep<1 or a.tol<0 or a.eta<=0: raise SystemExit('invalid parameters')
    if a.device=='mps' and not torch.backends.mps.is_available(): raise SystemExit('MPS unavailable')
    raw,candidate,C,factors,product_error=load_v2(a.input,a.candidate_index); g=len(factors); total=1<<(g-1); end=total if a.end is None else min(total,a.end)
    if not 0<=a.start<=end<=total: raise SystemExit('invalid half-open mask interval')
    if a.smoke and (a.device!='cpu' or end-a.start>1024 or a.batch_size>1024): raise SystemExit('--smoke requires CPU, <=1024 masks and batch')
    device=torch.device(a.device); dtype=torch.float64 if a.device=='cpu' else torch.float32; D=len(C)-1; start_time=perf_counter()
    factors_t=[torch.tensor(x,dtype=dtype,device=device) for x in factors]; degree_weight=torch.tensor([len(x)-1 for x in factors[:-1]],dtype=torch.int64)
    source_hash=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(); dependency_hash=hashlib.sha256(V1.read_bytes()).hexdigest()
    counts={'codes_seen':0,'constant_A_excluded':0,'normalization_failed':0,'fractional_masks':0,'closed_tolerance_fractional_masks':0,'ill_conditioned_fractional_masks':0}; best=[]; smoke=[]
    report={'format':'gpu-root-group-search-checkpoint-v2','status':'running','input_sha256':hashlib.sha256(a.input.read_bytes()).hexdigest(),'source_sha256':source_hash,'v1_helper_sha256':dependency_hash,
      'prepared_candidate':{'index':a.candidate_index,'name':candidate.get('name'),'degree':candidate.get('degree'),'group_count':g,'group_product_relative_error':product_error},
      'parameters':{'start':a.start,'end':end,'batch_size':a.batch_size,'device':a.device,'dtype':str(dtype),'tol':a.tol,'eta':a.eta,'condition_limit':a.condition_limit,'smoke':a.smoke,'group_limit_v2':30},
      'versions':{'torch':torch.__version__,'python':sys.version,'platform':platform.platform()},'counts':counts,'completed_mask_half_open':[a.start,a.start],'candidates':[],
      'claim_boundary':'Numerical ranking only: no float32/float64 non-hit is an exclusion; incomplete, near-boundary, and ill-conditioned masks are unknown.'}
    try:
      for lo in range(a.start,end,a.batch_size):
        hi=min(end,lo+a.batch_size); code_cpu=torch.arange(lo,hi,dtype=torch.int64); bits_cpu=((code_cpu[:,None]>>torch.arange(g-1))&1).bool(); da_cpu=(bits_cpu.to(torch.int64)*degree_weight).sum(1); db_cpu=D-da_cpu
        code=code_cpu.to(device); bits=bits_cpu.to(device); da=da_cpu.to(device); db=db_cpu.to(device); n=hi-lo
        A=torch.zeros((n,D+1),dtype=dtype,device=device); B=torch.zeros_like(A); A[:,0]=B[:,0]=1
        for j in range(g):
          select=bits[:,j] if j<g-1 else torch.zeros(n,dtype=torch.bool,device=device); MA,MB=tmul(A,factors_t[j]),tmul(B,factors_t[j]); A=torch.where(select[:,None],MA,A); B=torch.where(select[:,None],B,MB)
        nonconstant=code!=0; a0,b0=A[:,0],B[:,0]; normal=nonconstant&(a0.abs()>a.tol)&(b0.abs()>a.tol)
        counts['codes_seen']+=n; counts['constant_A_excluded']+=int((~nonconstant).sum().item()); counts['normalization_failed']+=int((nonconstant&~normal).sum().item())
        A=A/a0[:,None]; B=B/b0[:,None]; activeA=(torch.arange(D+1)[None,:]<=da_cpu[:,None]).to(device); activeB=(torch.arange(D+1)[None,:]<=db_cpu[:,None]).to(device); actual=torch.cat((torch.where(activeA,A,torch.zeros_like(A)),torch.where(activeB,B,torch.zeros_like(B))),1)
        mn=actual.min(1).values; mx=actual.max(1).values; range_violation=torch.relu(-actual-a.tol).sum(1)+torch.relu(actual-(1+a.tol)).sum(1); fractional=(((A>a.eta)&(A<1-a.eta)&activeA).any(1)|((B>a.eta)&(B<1-a.eta)&activeB).any(1)); maxabs=torch.maximum(torch.where(activeA,A.abs(),torch.zeros_like(A)).max(1).values,torch.where(activeB,B.abs(),torch.zeros_like(B)).max(1).values); cond=maxabs; leading=(A.gather(1,da[:,None]).squeeze(1)-1).abs()+(B.gather(1,db[:,None]).squeeze(1)-1).abs(); score=range_violation+leading; eligible=normal&fractional
        counts['fractional_masks']+=int(eligible.sum().item()); counts['closed_tolerance_fractional_masks']+=int((eligible&(mn>=-a.tol)&(mx<=1+a.tol)).sum().item()); counts['ill_conditioned_fractional_masks']+=int((eligible&(cond>a.condition_limit)).sum().item())
        values,pos=torch.topk(torch.where(eligible,score,torch.full_like(score,float('inf'))),k=min(a.keep,n),largest=False); keepmask=torch.isfinite(values); pos=pos[keepmask]
        # Exactly one device->CPU transfer for all selected numerical metrics.
        packed=torch.stack((score[pos],mn[pos],mx[pos],leading[pos],cond[pos]),1).detach().cpu().numpy(); chosen=pos.detach().cpu().numpy()
        for p,row in zip(chosen,packed): best.append({'mask':int(lo+int(p)),'float_device_score':float(row[0]),'float_device_min_actual':float(row[1]),'float_device_max_actual':float(row[2]),'float_device_leading_error':float(row[3]),'float_device_condition_proxy':float(row[4])})
        best=sorted({x['mask']:x for x in best}.values(),key=lambda x:(x['float_device_score'],x['mask']))[:a.keep]
        if a.smoke:
          for p in range(min(4,n)):
            if lo+p:
              NA,NB,_=v1.np_factors(factors,lo+p); smoke.append({'mask':lo+p,'torch_vs_numpy_max_abs':float(max(np.abs(A[p,:len(NA)].detach().cpu().numpy()-NA).max(),np.abs(B[p,:len(NB)].detach().cpu().numpy()-NB).max()))})
        del A,B,bits,code,actual,activeA,activeB
        report.update({'counts':counts,'completed_mask_half_open':[a.start,hi],'candidates':best,'smoke_samples':smoke[-32:],'elapsed_seconds':perf_counter()-start_time}); atomic(a.checkpoint,report)
      checked=[]
      for row in best: row['cpu_float64_recheck']=v1.recheck(C,factors,row['mask'],a.tol,a.eta,a.condition_limit); checked.append(row)
      report.update({'status':'complete','completed_mask_half_open':[a.start,end],'candidates':checked,'elapsed_seconds':perf_counter()-start_time}); atomic(a.checkpoint,report)
      print(json.dumps({'status':'complete','completed_mask_half_open':[a.start,end],'counts':counts,'retained':len(checked),'checkpoint':str(a.checkpoint)},indent=2))
    finally:
      if a.device=='mps': torch.mps.synchronize(); torch.mps.empty_cache()

if __name__=='__main__': main()
