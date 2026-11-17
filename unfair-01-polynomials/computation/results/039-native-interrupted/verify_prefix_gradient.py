#!/usr/bin/env python3
"""Float64 check of native-prefix-gradient.wl DFT/correlation backprop."""
from __future__ import annotations
import hashlib,json,sys
from pathlib import Path
import numpy as np
SEED=2026111617
def mats(q,l):
 r=np.arange(q)[:,None];s=np.arange(l)[None,:];t=2*np.pi*r*s/l
 return np.cos(t),-np.sin(t),np.cos(t).T/l,-np.sin(t).T/l
def sig(x):return 1/(1+np.exp(-x))
def calc(za,zb,m=6,n=9,l=32):
 faR,faI,iaR,iaI=mats(m+1,l);fbR,fbI,ibR,ibI=mats(n+1,l);fcR,fcI,icR,icI=mats(m+n+1,l)
 fa=np.zeros(m+1);fa[[0,1,m]]=[1,.5,1];fb=np.zeros(n+1);fb[[0,1,n]]=[1,.5,1];ma=1-(fa!=0);mb=1-(fb!=0)
 a=fa+ma*sig(za);b=fb+mb*sig(zb);ar=a@faR;ai=a@faI;br=b@fbR;bi=b@fbI;cr=ar*br-ai*bi;ci=ar*bi+ai*br;c=cr@icR+ci@icI
 e=2*c*(c-1)*(2*c-1);er=e@fcR;ei=e@fcI
 gra=(er*br+ei*bi)@iaR+(ei*br-er*bi)@iaI;grb=(er*ar+ei*ai)@ibR+(ei*ar-er*ai)@ibI
 loss=np.sum(c*c*(c-1)**2,axis=1);return loss,a,b,c,np.array([np.convolve(x,y) for x,y in zip(a,b)]),gra*ma*sig(za)*(1-sig(za)),grb*mb*sig(zb)*(1-sig(zb))
def main():
 rng=np.random.default_rng(SEED);za=rng.normal(size=(3,7));zb=rng.normal(size=(3,10));loss,a,b,c,direct,ga,gb=calc(za,zb)
 assert np.max(abs(c-direct))<1e-12 and np.all(ga[:,[0,1,6]]==0) and np.all(gb[:,[0,1,9]]==0)
 h=1e-6;checks=[]
 for w,ij in [('a',(1,3)),('a',(2,5)),('b',(0,4)),('b',(2,8))]:
  ap,bp,am,bm=za.copy(),zb.copy(),za.copy(),zb.copy()
  if w=='a':ap[ij]+=h;am[ij]-=h;g=ga[ij]
  else:bp[ij]+=h;bm[ij]-=h;g=gb[ij]
  fd=(calc(ap,bp)[0][ij[0]]-calc(am,bm)[0][ij[0]])/(2*h);checks.append({'coordinate':w+str(ij),'analytic':float(g),'finite_difference':float(fd),'absolute_error':float(abs(g-fd))})
 assert max(x['absolute_error'] for x in checks)<2e-7
 src=Path(__file__).with_name('native-prefix-gradient.wl');out={'status':'passed_float64_dft_correlation_gradient_check','seed':SEED,'small_model':{'m':6,'n':9,'L':32,'starts':3},'max_linear_convolution_error':float(np.max(abs(c-direct))),'fixed_logit_gradients_zero':True,'finite_difference_checks':checks,'max_finite_difference_error':max(x['absolute_error'] for x in checks),'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'python':sys.version,'numpy':np.__version__,'scope':'Formula check only; no Wolfram GPU dispatch/global optimization claim.'}
 Path(__file__).with_name('verify_prefix_gradient.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps(out,sort_keys=True))
if __name__=='__main__':main()
