#!/usr/bin/env python3
"""Direct polynomial CPU cross-check for the native six-mask Metal smoke."""
import hashlib,json,sys
from pathlib import Path
import numpy as np
ROOT=Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/round038-macbook-gpu52-direct')
def q(g):
 z=np.array([complex(*x) for x in g]);p=np.polynomial.polynomial.polyfromroots(z);p=p/p[0];assert np.max(abs(p.imag))<1e-7;return p.real
def pair(t,m):
 bits=[0]+list(map(int,format(m,f'0{t["groups"]-1}b')));a=np.array([1.]);b=np.array([1.])
 for bit,g in zip(bits,t['roots']):
  if bit:a=np.polynomial.polynomial.polymul(a,q(g))
  else:b=np.polynomial.polynomial.polymul(b,q(g))
 return np.pad(a,(0,53-len(a))),np.pad(b,(0,53-len(b)))
def score(a,b):return max(0.,-a.min(),-b.min(),a.max()-1,b.max()-1)
out=Path(sys.argv[1]);got=json.loads(out.read_text());d=json.loads((ROOT/'inputs.json').read_text());t=d['tasks'][2];err=[]
for i,m in enumerate(got['masks']):
 a,b=pair(t,m);want=score(a,b);g=got['scores'][i]
 if m==0: assert g is None;continue
 g=float(g)
 err.append(abs(g-want))
assert all(x==1 for x in got['finite']) and max(err)<1e-3,(max(err),got)
report={'status':'passed_native_metal_six_mask_cpu_polynomial_crosscheck','scope':'Float32 GPU ranking agreement only; CPU direct polynomial reconstruction remains the candidate verifier, not a certificate.','masks':got['masks'],'max_score_error':max(err),'metal_device':got['metal_device'],'dims':{'groups':got['groups'],'nfreq':got['nfreq'],'ncoef':got['ncoef']},'top_indices':got['top_indices'],'top_scores':got['top_scores'],'source_sha256':hashlib.sha256(Path(__file__).with_name('NativeScore.swift').read_bytes()).hexdigest(),'output_sha256':hashlib.sha256(out.read_bytes()).hexdigest()}
Path(__file__).with_name('cpu-verify.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');print(json.dumps(report,sort_keys=True))
