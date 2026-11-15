#!/usr/bin/env python3
"""Tiny actual-MPS fair-mask calibration for double-single root products."""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import torch

SRC=r'''#pragma clang fp contract(off)
#include <metal_stdlib>
using namespace metal;
inline float2 qsum(float a,float b){float s=a+b;return float2(s,b-(s-a));}
inline float2 tsum(float a,float b){float s=a+b,bb=s-a;return float2(s,(a-(s-bb))+(b-bb));}
inline float2 tprod(float a,float b){float p=a*b;return float2(p,fma(a,b,-p));}
inline float2 dsadd(float2 a,float2 b){float2 s=tsum(a.x,b.x);return qsum(s.x,s.y+a.y+b.y);}
inline float2 dsmul(float2 a,float2 b){float2 p=tprod(a.x,b.x);return qsum(p.x,p.y+a.x*b.y+a.y*b.x+a.y*b.y);}
kernel void dsconv(device float* oh [[buffer(0)]],device float* ol [[buffer(1)]],device const float* ah [[buffer(2)]],device const float* al [[buffer(3)]],device const float* bh [[buffer(4)]],device const float* bl [[buffer(5)]],constant uint& aw [[buffer(6)]],constant uint& bw [[buffer(7)]],uint k [[thread_position_in_grid]]){
 if(k>=aw+bw-1)return;uint lo=(k>=bw-1)?k-(bw-1):0,hi=min(k,aw-1);float2 z=float2(0.0f);for(uint i=lo;i<=hi;i++)z=dsadd(z,dsmul(float2(ah[i],al[i]),float2(bh[k-i],bl[k-i])));oh[k]=z.x;ol[k]=z.y;
}'''
def split(x):
 h=np.asarray(x,np.float32);return h,np.asarray(np.asarray(x,np.float64)-h.astype(np.float64),np.float32)
def run_product(lib,groups,indices):
 h,l=np.array([1.],np.float32),np.array([0.],np.float32)
 for j in indices:
  bh,bl=split(groups[j]);aw,bw=len(h),len(bh);out=aw+bw-1
  # First tensor is output: compile_shader's inferred grid has exactly out threads.
  oh=torch.empty(out,device='mps');ol=torch.empty(out,device='mps');ah=torch.from_numpy(np.pad(h,(0,out-aw))).to('mps');al=torch.from_numpy(np.pad(l,(0,out-aw))).to('mps');tb=torch.from_numpy(bh).to('mps');tl=torch.from_numpy(bl).to('mps')
  lib.dsconv(oh,ol,ah,al,tb,tl,aw,bw);torch.mps.synchronize();h=oh.cpu().numpy();l=ol.cpu().numpy()
 return h.astype(np.float64)+l.astype(np.float64)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('input',type=Path);ap.add_argument('--candidate-index',type=int,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
 if not torch.backends.mps.is_available():raise SystemExit('MPS unavailable')
 data=json.loads(a.input.read_text());c=data['candidates'][a.candidate_index];L=c['L'];groups=[np.asarray(x['coefficients_ascending'],float)for x in c['groups']];original={int(b['original_group'])for b in c['blocks']};lifts=set(range(len(groups)))-original;H=np.asarray(c['H_coefficients_ascending'],float);HL=np.zeros((len(H)-1)*L+1);HL[::L]=H
 try:
  lib=torch.mps.compile_shader(SRC)
  A=run_product(lib,groups,sorted(original));B=run_product(lib,groups,sorted(lifts))
  # The balanced no-real-root identity makes true endpoints one.  Set these
  # exact endpoints before evaluating interior cancellation error.
  A[0]=A[-1]=1.;B[0]=B[-1]=1.
  out={'format':'mps-double-single-fair-calibration-v1','candidate':c['name'],'L':L,'device':'mps','torch':torch.__version__,'original_group_indices':sorted(original),'lift_group_count':len(lifts),'H_max_abs_error':float(abs(A-H).max()),'HxL_max_abs_error':float(abs(B-HL).max()),'fair_accept_lt_1e-9':bool(max(abs(A-H).max(),abs(B-HL).max())<1e-9),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'input_sha256':hashlib.sha256(a.input.read_bytes()).hexdigest(),'metal_contract_request':'#pragma clang fp contract(off)','claim_boundary':'Actual MPS numerical calibration only; not a proof or an exclusion.'}
  a.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
 finally:
  torch.mps.synchronize();torch.mps.empty_cache()
if __name__=='__main__':main()
