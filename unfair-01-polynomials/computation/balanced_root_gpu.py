#!/usr/bin/env python3
"""Balanced root-group numerical screen; CPU default, MPS explicitly opt-in.

Expected prepared JSON has L, C_coefficients_ascending (or
coefficients_ascending), groups_ascending_coefficients, and blocks.  A block
has kind, original_group (or original_group_index), and lift_groups (or
lift_group_indices).  Group indices refer to the ascending real factors.
"""
from __future__ import annotations
import argparse, hashlib, itertools, json, platform, sys
from pathlib import Path
from time import perf_counter
import numpy as np
import torch

def atomic(p,x):
 t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(x,indent=2)+'\n');t.replace(p)
def tmul(p,f):
 o=torch.zeros_like(p);w=p.shape[1]
 for s in range(f.numel()):o[:,s:]+=f[s]*p[:,:w-s]
 return o
def fields(b): return int(b.get('original_group',b.get('original_group_index'))),[int(x) for x in b.get('lift_groups',b.get('lift_group_indices',[]))]
def local_options(block,L,E,weights):
 original,lifts=fields(block)
 if E==L:return [1<<original,(1<<original)*0+sum(1<<j for j in lifts)]
 out=[]
 for bits in range(1<<len(lifts)):
  if sum(weights[lifts[q]] for q in range(len(lifts)) if(bits>>q)&1)==E:out.append(sum(1<<lifts[q] for q in range(len(lifts)) if(bits>>q)&1))
 return out
def npcheck(C,fs,mask,tol,eta,clim):
 A=np.array([1.]);B=np.array([1.]);da=0
 for j,f in enumerate(fs):
  if(mask>>j)&1:A=np.convolve(A,f);da+=len(f)-1
  else:B=np.convolve(B,f)
 A=A/A[0];B=B/B[0];z=np.r_[A,B];conv=np.convolve(A,B);return {'mask':int(mask),'degree_A':da,'degree_B':len(C)-1-da,'A_coefficients_ascending':A.tolist(),'B_coefficients_ascending':B.tolist(),'convolution_max_abs_residual':float(abs(conv-C/C[0]).max()),'minimum_actual_coefficient':float(z.min()),'maximum_actual_coefficient':float(z.max()),'within_closed_tolerance':bool(z.min()>=-tol and z.max()<=1+tol),'has_fractionality':bool(((z>eta)&(z<1-eta)).any()),'near_boundary':bool((abs(z)<=10*tol).any()or(abs(z-1)<=10*tol).any()),'condition_proxy':float(abs(z).max()),'ill_conditioned':bool(abs(z).max()>clim),'normalized_leading_error':float(abs(A[-1]-1)+abs(B[-1]-1))}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('input',type=Path);ap.add_argument('--candidate-index',type=int,default=0);ap.add_argument('--checkpoint',type=Path,required=True);ap.add_argument('--device',choices=('cpu','mps'),default='cpu');ap.add_argument('--batch-size',type=int,default=131072);ap.add_argument('--max-codes',type=int);ap.add_argument('--tol',type=float,default=1e-6);ap.add_argument('--eta',type=float,default=1e-5);ap.add_argument('--condition-limit',type=float,default=1e8);ap.add_argument('--smoke',action='store_true');a=ap.parse_args()
 container=json.loads(a.input.read_text());candidates=container.get('candidates');raw=(candidates[a.candidate_index] if isinstance(candidates,list) else container)
 if not 0<=a.candidate_index<(len(candidates) if isinstance(candidates,list) else 1):raise SystemExit('invalid candidate index')
 L=int(raw['L']);C=np.asarray(raw.get('C_coefficients_ascending',raw.get('coefficients_ascending')),dtype=np.float64);groups_raw=raw.get('groups');fs=[np.asarray(q['coefficients_ascending'],dtype=np.float64)for q in groups_raw] if groups_raw is not None else [np.asarray(q,dtype=np.float64)for q in raw['groups_ascending_coefficients']];blocks=raw['blocks'];g=len(fs)
 if not(2<=g<=48 and L>=2 and C.ndim==1 and abs(C[0])>1e-14):raise SystemExit('invalid L/groups/C')
 if a.device=='mps'and not torch.backends.mps.is_available():raise SystemExit('MPS unavailable')
 if a.smoke and(a.device!='cpu'or a.max_codes is None or a.max_codes>1024 or a.batch_size>1024):raise SystemExit('--smoke requires CPU and max-codes,batch <=1024')
 product=np.array([1.]);
 for f in fs:
  if len(f)not in(2,3)or f[-1]!=1:raise SystemExit('groups must be monic real linear/quadratic')
  product=np.convolve(product,f)
 perr=float(abs(product-C).max()/max(1.,abs(C).max()))
 if perr>1e-8:raise SystemExit('group product does not reproduce C')
 bw=np.asarray([q['weight'] for q in groups_raw],dtype=np.int64) if groups_raw is not None else np.asarray([len(f)-1 for f in fs],dtype=np.int64)
 if bw.shape!=(g,)or(bw<=0).any():raise SystemExit('invalid group balance weights')
 device=torch.device(a.device);dtype=torch.float64 if a.device=='cpu'else torch.float32;D=len(C)-1;ft=[torch.tensor(f,dtype=dtype,device=device)for f in fs];dw=np.asarray([len(f)-1 for f in fs],dtype=np.int64);src=hashlib.sha256(Path(__file__).read_bytes()).hexdigest();inp=hashlib.sha256(a.input.read_bytes()).hexdigest();started=perf_counter();progress=[]
 report={'format':'balanced-root-gpu-checkpoint-v1','status':'running','source_sha256':src,'input_sha256':inp,'prepared_candidate':{'index':a.candidate_index,'name':raw.get('name'),'total_groups':g},'parameters':{'L':L,'device':a.device,'dtype':str(dtype),'batch_size':a.batch_size,'max_codes':a.max_codes,'tol':a.tol,'eta':a.eta,'condition_limit':a.condition_limit,'group_limit':48},'group_count':g,'maximum_degree':D,'balance_weights':bw.tolist(),'group_product_relative_error':perr,'versions':{'torch':torch.__version__,'python':sys.version,'platform':platform.platform()},'E_progress':progress,'claim_boundary':'Numerical ranking only: incomplete E intervals, near-boundary, ill-conditioned, and every float non-hit are unknown.'}
 try:
  for E in range(1,L+1):
   opts=[local_options(b,L,E,bw)for b in blocks]
   if E==L:opts[0]=[1<<fields(blocks[0])[0]]
   if any(not x for x in opts):progress.append({'E':E,'status':'no_local_options'});atomic(a.checkpoint,report);continue
   rad=[len(x)for x in opts];total=1
   for r in rad:total*=r
   if total>(1<<63)-1:progress.append({'E':E,'status':'count_exceeds_int64','total_codes_decimal':str(total)});atomic(a.checkpoint,report);continue
   end=total if a.max_codes is None else min(total,a.max_codes);best=[];counts={'codes_seen':0,'fractional_masks':0,'closed_tolerance_fractional_masks':0,'ill_conditioned_fractional_masks':0};smoke=[]
   for lo in range(0,end,a.batch_size):
    hi=min(end,lo+a.batch_size);q=np.arange(lo,hi,dtype=np.int64);m=np.zeros(hi-lo,dtype=np.uint64)
    for options,r in zip(opts,rad):idx=q%r;q//=r;m|=np.asarray(options,dtype=np.uint64)[idx]
    bits_cpu=((m[:,None]>>np.arange(g,dtype=np.uint64))&1).astype(np.bool_);da_cpu=bits_cpu.astype(np.int64)@dw;db_cpu=D-da_cpu;n=hi-lo;bits=torch.from_numpy(bits_cpu).to(device);da=torch.from_numpy(da_cpu).to(device);db=torch.from_numpy(db_cpu).to(device)
    A=torch.zeros((n,D+1),dtype=dtype,device=device);B=torch.zeros_like(A);A[:,0]=B[:,0]=1
    for j in range(g):sel=bits[:,j];MA,MB=tmul(A,ft[j]),tmul(B,ft[j]);A=torch.where(sel[:,None],MA,A);B=torch.where(sel[:,None],B,MB)
    a0,b0=A[:,0],B[:,0];normal=(a0.abs()>a.tol)&(b0.abs()>a.tol);A=A/a0[:,None];B=B/b0[:,None];activeA=(torch.arange(D+1)[None,:]<=torch.from_numpy(da_cpu)[:,None]).to(device);activeB=(torch.arange(D+1)[None,:]<=torch.from_numpy(db_cpu)[:,None]).to(device);actual=torch.cat((torch.where(activeA,A,torch.zeros_like(A)),torch.where(activeB,B,torch.zeros_like(B))),1);mn=actual.min(1).values;mx=actual.max(1).values;rv=torch.relu(-actual-a.tol).sum(1)+torch.relu(actual-(1+a.tol)).sum(1);frac=(((A>a.eta)&(A<1-a.eta)&activeA).any(1)|((B>a.eta)&(B<1-a.eta)&activeB).any(1));maxabs=torch.maximum(torch.where(activeA,A.abs(),torch.zeros_like(A)).max(1).values,torch.where(activeB,B.abs(),torch.zeros_like(B)).max(1).values);lead=(A.gather(1,da[:,None]).squeeze(1)-1).abs()+(B.gather(1,db[:,None]).squeeze(1)-1).abs();score=rv+lead;ok=normal&frac;counts['codes_seen']+=n;counts['fractional_masks']+=int(ok.sum().item());counts['closed_tolerance_fractional_masks']+=int((ok&(mn>=-a.tol)&(mx<=1+a.tol)).sum().item());counts['ill_conditioned_fractional_masks']+=int((ok&(maxabs>a.condition_limit)).sum().item());val,pos=torch.topk(torch.where(ok,score,torch.full_like(score,float('inf'))),k=min(64,n),largest=False);pos=pos[torch.isfinite(val)];pack=torch.stack((score[pos],mn[pos],mx[pos],lead[pos],maxabs[pos]),1).detach().cpu().numpy();pc=pos.detach().cpu().numpy()
    for p,row in zip(pc,pack):best.append({'mask':int(m[p]),'code':int(lo+p),'float_device_score':float(row[0]),'float_device_min_actual':float(row[1]),'float_device_max_actual':float(row[2]),'float_device_leading_error':float(row[3]),'float_device_condition_proxy':float(row[4])})
    best=sorted({x['mask']:x for x in best}.values(),key=lambda x:(x['float_device_score'],x['mask']))[:64]
    if a.smoke:
     for p in range(min(4,n)):
      checked=npcheck(C,fs,int(m[p]),a.tol,a.eta,a.condition_limit);NA=np.asarray(checked['A_coefficients_ascending']);NB=np.asarray(checked['B_coefficients_ascending']);err=max(abs(A[p,:len(NA)].detach().cpu().numpy()-NA).max(),abs(B[p,:len(NB)].detach().cpu().numpy()-NB).max());smoke.append({'mask':int(m[p]),'torch_vs_numpy_max_abs':float(err)})
    del A,B,bits,actual
    entry={'E':E,'status':'running','total_codes':total,'completed_code_half_open':[0,hi],'counts':counts,'top64':best,'smoke_samples':smoke[-32:]};old=[x for x in progress if x.get('E')!=E];progress[:]=old+[entry];report['elapsed_seconds']=perf_counter()-started;atomic(a.checkpoint,report)
   for z in best:z['cpu_float64_recheck']=npcheck(C,fs,z['mask'],a.tol,a.eta,a.condition_limit)
   entry={'E':E,'status':'complete'if end==total else'partial_unknown','total_codes':total,'completed_code_half_open':[0,end],'counts':counts,'top64':best,'smoke_samples':smoke[-32:]};progress[:]=[x for x in progress if x.get('E')!=E]+[entry];atomic(a.checkpoint,report)
  report['status']='complete'if all(x.get('status')in('complete','no_local_options')for x in progress)else'partial_unknown';report['elapsed_seconds']=perf_counter()-started;atomic(a.checkpoint,report);print(json.dumps({'status':report['status'],'E_count':len(progress),'checkpoint':str(a.checkpoint)},indent=2))
 finally:
  if a.device=='mps':torch.mps.synchronize();torch.mps.empty_cache()
if __name__=='__main__':main()
