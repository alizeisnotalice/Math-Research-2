#!/usr/bin/env python3
"""Checkpointed PyTorch root-group mask screen (default CPU; MPS opt-in).

Consumes algebraic_global/root_group_candidates_20261115.json.  Each candidate
has ascending C coefficients and real group factors, of length two (linear)
or three (a conjugate-pair quadratic).  The final group is permanently B;
bits 0..g-2 select groups for A.  This is numerical reconnaissance only.
"""
from __future__ import annotations
import argparse, hashlib, json, platform, sys
from pathlib import Path
from time import perf_counter
from typing import Any
import numpy as np
import torch


def die(message: str) -> None: raise SystemExit(message)

def atomic(path: Path, item: dict[str, Any]) -> None:
    tmp = path.with_suffix(path.suffix+'.tmp'); tmp.write_text(json.dumps(item, indent=2)+'\n'); tmp.replace(path)

def load(path: Path, candidate_index: int):
    raw = json.loads(path.read_text()); candidates = raw.get('candidates')
    if not isinstance(candidates, list) or not 0 <= candidate_index < len(candidates): die('invalid candidate index')
    c = candidates[candidate_index]
    C = np.asarray(c.get('coefficients_ascending'), dtype=np.float64)
    factors = c.get('groups_ascending_coefficients')
    if C.ndim != 1 or len(C) < 2 or not isinstance(factors, list) or not 2 <= len(factors) <= 25: die('invalid prepared candidate')
    fs = [np.asarray(x, dtype=np.float64) for x in factors]
    if any(f.ndim != 1 or len(f) not in (2, 3) or f[-1] != 1.0 for f in fs): die('groups must be monic real degree-1/2 ascending factors')
    product = np.array([1.0])
    for f in fs: product = np.convolve(product, f)
    err = np.max(np.abs(product-C))/max(1.0, np.max(np.abs(C)))
    if err > 1e-8 or abs(C[0]) < 1e-14: die(f'prepared product mismatch ({err}) or zero C(0)')
    return raw, c, C, fs, err

def tmul(poly: torch.Tensor, factor: torch.Tensor) -> torch.Tensor:
    out = torch.zeros_like(poly); width = poly.shape[1]
    for shift in range(factor.numel()): out[:, shift:] += factor[shift]*poly[:, :width-shift]
    return out

def np_factors(factors: list[np.ndarray], mask: int):
    A = np.array([1.0]); B = np.array([1.0]); da = 0
    for j, f in enumerate(factors):
        if j < len(factors)-1 and ((mask >> j) & 1): A=np.convolve(A,f); da += len(f)-1
        else: B=np.convolve(B,f)
    return A/A[0], B/B[0], da

def recheck(C, factors, mask, tol, eta, condition_limit):
    A,B,da=np_factors(factors,mask); DB=len(C)-1-da; conv=np.convolve(A,B); target=C/C[0]
    z=np.r_[A,B]; mn=float(z.min()); mx=float(z.max()); cond=float(max(abs(A).max(),abs(B).max())/min(abs(A[0]),abs(B[0])))
    leading_a,leading_b=float(A[-1]),float(B[-1])
    return {'mask':int(mask),'degree_A':da,'degree_B':DB,'A_coefficients_ascending':A.tolist(),'B_coefficients_ascending':B.tolist(),
            'convolution_max_abs_residual':float(abs(conv-target).max()),'minimum_actual_coefficient':mn,'maximum_actual_coefficient':mx,
            'within_closed_tolerance':bool(mn>=-tol and mx<=1+tol),'has_fractionality':bool(((z>eta)&(z<1-eta)).any()),
            'normalized_leading_A':leading_a,'normalized_leading_B':leading_b,
            'normalized_leading_error':abs(leading_a-1)+abs(leading_b-1),
            'leading_within_tolerance':bool(abs(leading_a-1)<=tol and abs(leading_b-1)<=tol),
            'near_boundary':bool((abs(z)<=10*tol).any() or (abs(z-1)<=10*tol).any()),'condition_proxy':cond,'ill_conditioned':bool(cond>condition_limit)}

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('input',type=Path); ap.add_argument('--candidate-index',type=int,required=True)
    ap.add_argument('--checkpoint',type=Path,required=True); ap.add_argument('--start',type=int,default=0); ap.add_argument('--end',type=int)
    ap.add_argument('--batch-size',type=int,default=65536); ap.add_argument('--keep',type=int,default=64); ap.add_argument('--device',choices=('cpu','mps'),default='cpu')
    ap.add_argument('--tol',type=float,default=1e-6); ap.add_argument('--eta',type=float,default=1e-5); ap.add_argument('--condition-limit',type=float,default=1e8); ap.add_argument('--smoke',action='store_true')
    a=ap.parse_args()
    if a.batch_size<1 or a.keep<1 or a.tol<0 or a.eta<=0: die('invalid numerical parameters')
    if a.device=='mps' and not torch.backends.mps.is_available(): die('MPS unavailable')
    if a.smoke and (a.device!='cpu' or (a.end is not None and a.end-a.start>1024) or a.batch_size>1024): die('--smoke requires CPU and an interval/batch no larger than 1024')
    raw,candidate,C,factors,product_error=load(a.input,a.candidate_index); g=len(factors); allcodes=1<<(g-1); end=allcodes if a.end is None else min(a.end,allcodes)
    if not 0<=a.start<=end<=allcodes: die('invalid half-open mask interval')
    if a.smoke and end-a.start>1024: die('--smoke interval is larger than 1024 masks')
    device=torch.device(a.device); dtype=torch.float64 if a.device=='cpu' else torch.float32; D=len(C)-1; started=perf_counter()
    f=[torch.tensor(x,dtype=dtype,device=device) for x in factors]; degree_weights=torch.tensor([len(x)-1 for x in factors[:-1]],dtype=torch.int64)
    counts={'codes_seen':0,'constant_A_excluded':0,'normalization_failed':0,'fractional_masks':0,'closed_tolerance_fractional_masks':0,'ill_conditioned_fractional_masks':0}; best=[]; smoke=[]
    report={'format':'gpu-root-group-search-checkpoint-v1','status':'running','input_sha256':hashlib.sha256(a.input.read_bytes()).hexdigest(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'prepared_candidate':{'index':a.candidate_index,'name':candidate.get('name'),'degree':candidate.get('degree'),'group_count':g,'group_product_relative_error':product_error},
            'parameters':{'start':a.start,'end':end,'batch_size':a.batch_size,'device':a.device,'dtype':str(dtype),'tol':a.tol,'eta':a.eta,'condition_limit':a.condition_limit,'smoke':a.smoke},
            'versions':{'torch':torch.__version__,'python':sys.version,'platform':platform.platform()},'counts':counts,'completed_mask_half_open':[a.start,a.start],'candidates':[],
            'claim_boundary':'A float32/float64 scan is not an exclusion theorem. Unscanned, conditioning-flagged, and near-boundary masks remain unknown.'}
    try:
      for lo in range(a.start,end,a.batch_size):
        hi=min(lo+a.batch_size,end)
        # Bit arithmetic remains on CPU: MPS support for int64 shifts varies
        # by PyTorch release.  Only the polynomial recurrence is device work.
        code_cpu=torch.arange(lo,hi,dtype=torch.int64); bits_cpu=((code_cpu[:,None]>>torch.arange(g-1))&1).bool(); da_cpu=(bits_cpu.to(torch.int64)*degree_weights).sum(1); db_cpu=D-da_cpu
        code=code_cpu.to(device); bits=bits_cpu.to(device); da=da_cpu.to(device); db=db_cpu.to(device)
        A=torch.zeros((hi-lo,D+1),dtype=dtype,device=device); B=torch.zeros_like(A); A[:,0]=B[:,0]=1
        for j in range(g):
          take=bits[:,j] if j<g-1 else torch.zeros(hi-lo,dtype=torch.bool,device=device)
          MA,MB=tmul(A,f[j]),tmul(B,f[j]); A=torch.where(take[:,None],MA,A); B=torch.where(take[:,None],B,MB)
        nonconstant=code!=0; a0,b0=A[:,0],B[:,0]; normal=nonconstant&(a0.abs()>a.tol)&(b0.abs()>a.tol); counts['codes_seen']+=hi-lo; counts['constant_A_excluded']+=int((~nonconstant).sum().item()); counts['normalization_failed']+=int((nonconstant&~normal).sum().item())
        A=A/a0[:,None]; B=B/b0[:,None]
        activeA=(torch.arange(D+1)[None,:]<=da_cpu[:,None]).to(device); activeB=(torch.arange(D+1)[None,:]<=db_cpu[:,None]).to(device); actual=torch.cat((torch.where(activeA,A,torch.zeros_like(A)),torch.where(activeB,B,torch.zeros_like(B))),1)
        mn=actual.min(1).values; mx=actual.max(1).values; range_violation=torch.relu(-actual-a.tol).sum(1)+torch.relu(actual-(1+a.tol)).sum(1); fractional=(((A>a.eta)&(A<1-a.eta)&activeA).any(1)|((B>a.eta)&(B<1-a.eta)&activeB).any(1)); cond=torch.maximum(torch.where(activeA,A.abs(),torch.zeros_like(A)).max(1).values,torch.where(activeB,B.abs(),torch.zeros_like(B)).max(1).values)
        leadA=A.gather(1,da[:,None]).squeeze(1); leadB=B.gather(1,db[:,None]).squeeze(1); leading_error=(leadA-1).abs()+(leadB-1).abs(); score=range_violation+leading_error; good=normal&fractional
        counts['fractional_masks']+=int(good.sum().item()); counts['closed_tolerance_fractional_masks']+=int((good&(mn>=-a.tol)&(mx<=1+a.tol)).sum().item()); counts['ill_conditioned_fractional_masks']+=int((good&(cond>a.condition_limit)).sum().item())
        vals,pos=torch.topk(torch.where(good,score,torch.full_like(score,float('inf'))),k=min(a.keep,hi-lo),largest=False)
        for v,p in zip(vals.detach().cpu().tolist(),pos.detach().cpu().tolist()):
          if np.isfinite(v): best.append({'mask':int(lo+p),'float_device_score':float(v),'float_device_min_actual':float(mn[p].cpu()),'float_device_max_actual':float(mx[p].cpu()),'float_device_leading_error':float(leading_error[p].cpu()),'float_device_condition_proxy':float(cond[p].cpu())})
        best=sorted({x['mask']:x for x in best}.values(),key=lambda x:(x['float_device_score'],x['mask']))[:a.keep]
        if a.smoke:
          for p in range(min(4,hi-lo)):
            mask=int(lo+p)
            if mask:
              NA,NB,_=np_factors(factors,mask); ea=float(max(abs(A[p,:len(NA)].detach().cpu().numpy()-NA).max(),abs(B[p,:len(NB)].detach().cpu().numpy()-NB).max())); smoke.append({'mask':mask,'torch_vs_numpy_max_abs':ea})
        del A,B,bits,code,actual
        if a.device=='mps': torch.mps.synchronize(); torch.mps.empty_cache()
        report.update({'counts':counts,'completed_mask_half_open':[a.start,hi],'candidates':best,'smoke_samples':smoke[-32:],'elapsed_seconds':perf_counter()-started}); atomic(a.checkpoint,report)
      checked=[]
      for row in best:
        row['cpu_float64_recheck']=recheck(C,factors,row['mask'],a.tol,a.eta,a.condition_limit); checked.append(row)
      report.update({'status':'complete','completed_mask_half_open':[a.start,end],'candidates':checked,'elapsed_seconds':perf_counter()-started}); atomic(a.checkpoint,report)
      print(json.dumps({'status':'complete','completed_mask_half_open':[a.start,end],'counts':counts,'retained':len(checked),'checkpoint':str(a.checkpoint)},indent=2))
    finally:
      if a.device=='mps': torch.mps.synchronize(); torch.mps.empty_cache()

if __name__=='__main__': main()
