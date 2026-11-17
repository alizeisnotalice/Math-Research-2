#!/usr/bin/env python3
"""Exact weight histogram for all k=5 spread-seed transvections; no kernel scan."""
import hashlib,json,platform,sys,time
from pathlib import Path
K=5;N=1024;MOD=0b100101
def dot(a,b):return (a&b).bit_count()&1
def mul(a,b):
 o=0
 while b:
  if b&1:o^=a
  b>>=1;a<<=1
  if a&32:a^=MOD
 return o&31
def fwht(a):
 a=list(a);h=1
 while h<N:
  for i in range(0,N,2*h):
   for j in range(i,i+h):x,y=a[j],a[j+h];a[j],a[j+h]=x+y,x-y
  h*=2
 return a
def main():
 t=time.monotonic();neg=list(range(16));F=[]
 for z in range(N):
  x=z&31;y=z>>5;value=0
  for m in neg:
   if y==mul(m,x):value^=(z!=0)
  F.append(value)
 # Verify the base only once by its exact Walsh spectrum.
 assert set(abs(x) for x in fwht([1-2*b for b in F]))=={32}
 candidates=[];tested=0;hist={};walsh_hist={}
 for v in range(1,N):
  D=[F[x]^F[x^v] for x in range(N)];W=fwht([1-2*b for b in D]);S=W[0]
  # Every nonzero derivative of a bent function is balanced.
  assert S==0
  for lam in range(1,N):
   if dot(lam,v):continue
   tested+=1;weight=(N-S+W[lam])//4
   assert 4*weight==N-S+W[lam]
   hist[weight]=hist.get(weight,0)+1
   walsh_hist[W[lam]]=walsh_hist.get(W[lam],0)+1
   if weight==64:candidates.append({'v':v,'lambda':lam,'derivative_weight':(N-S)//2,'walsh_derivative_at_lambda':W[lam]})
 assert tested==522753
 assert sum(hist.values())==tested
 seed={'schema':'a28-r032-k5-spread-input-v1','k':5,'modulus_bin':bin(MOD),'negative_slopes':neg,'field_order':32}
 ip=Path(__file__).with_name('k5-transvection-input.json');ip.write_text(json.dumps(seed,indent=2)+'\n')
 out={'schema':'a28-r032-k5-transvection-histogram-v2','scope':'one R23 Desarguesian k=5 seed; every nonidentity invertible transvection T(x)=x+lambda(x)v, lambda(v)=0, lambda nonzero. This is a finite-family result only; no projection scan.','arithmetic':'exact GF(2) truth tables and integer FWHT','weight_formula':'|{x: lambda(x)=1, D_v F(x)=1}|=(N-W_{D_vF}(0)+W_{D_vF}(lambda))/4; W_{D_vF}(0)=0 here','input_path':str(ip),'input_sha256':hashlib.sha256(ip.read_bytes()).hexdigest(),'transvections_tested':tested,'weight_histogram':{str(k):hist[k] for k in sorted(hist)},'derivative_walsh_histogram_nonzero_lambda':{str(k):walsh_hist[k] for k in sorted(walsh_hist)},'minimum_positive_weight':min(k for k in hist if k>0),'maximum_derivative_walsh_nonzero_lambda':max(walsh_hist),'minimum_derivative_walsh_nonzero_lambda':min(walsh_hist),'weight64_candidate_count':len(candidates),'candidates':candidates,'elapsed_seconds':time.monotonic()-t,'python':sys.version,'platform':platform.platform(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 Path(__file__).with_name('k5-transvection-candidates.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k] for k in ['transvections_tested','weight_histogram','minimum_positive_weight','weight64_candidate_count','elapsed_seconds']},indent=2))
if __name__=='__main__':main()
