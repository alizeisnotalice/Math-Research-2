#!/usr/bin/env python3
"""L-balanced high-precision local block preparation for DS GPU combination.

Precomputes each parent block's (Ablock,Bblock) options at mp.dps=100.  A
future GPU pass need only choose one option per each of three blocks and make
two DS convolutions for A and two for B; it must not multiply 39 root groups.
"""
from __future__ import annotations
import argparse, hashlib, itertools, json
from pathlib import Path
import mpmath as mp
import numpy as np

def conv(a,b):
 o=[mp.mpf('0')]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):o[i+j]+=x*y
 return o
def pair(poly):
 x=np.asarray([float(v) for v in poly],np.float64);h=x.astype(np.float32);return {'hi':h.tolist(),'lo':(x-h.astype(np.float64)).astype(np.float32).tolist(),'float64':x.tolist()}
def quad(z):return [abs(z)**2,-2*mp.re(z),mp.mpf(1)]
def main():
 ap=argparse.ArgumentParser();ap.add_argument('input',type=Path);ap.add_argument('--candidate-index',type=int,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();mp.mp.dps=100
 raw=json.loads(a.input.read_text());c=raw['candidates'][a.candidate_index];L=int(c['L']);groups=c['groups'];blocks=c['blocks'];prepared=[]
 for block in blocks:
  re,im=map(mp.mpf,block['parent_root_positive_imag']);alpha=mp.mpc(re,im);r=abs(alpha)**(mp.mpf(1)/L);theta=mp.arg(alpha)
  generated=[quad(r*mp.e**(1j*(theta+2*mp.pi*k)/L)) for k in range(L)]
  lift_idx=list(map(int,block['lift_groups']));target=[np.asarray(groups[j]['coefficients_ascending'],float) for j in lift_idx]
  # Match high-precision lifts to the preparer's stable floating group order.
  unused=set(range(L));ordered=[]
  for t in target:
   k=min(unused,key=lambda q:float(np.max(np.abs(np.asarray([float(x) for x in generated[q]])-t))));unused.remove(k);ordered.append(generated[k])
  original=quad(alpha);opts={}
  for E in range(1,L):
   rows=[]
   for chosen in sorted(itertools.combinations(range(L),E),key=lambda s:sum(1<<q for q in s)):
    A=[mp.mpf(1)];B=original[:];chosen=set(chosen)
    for q,f in enumerate(ordered):
     if q in chosen:A=conv(A,f)
     else:B=conv(B,f)
    rows.append({'selected_lift_positions':sorted(chosen),'Ablock':pair(A),'Bblock':pair(B)})
   opts[str(E)]=rows
  # Product of all lift quadratics is exactly x^(2L)-2Re(alpha)x^L+|alpha|^2.
  full=[abs(alpha)**2]+[mp.mpf(0)]*(L-1)+[-2*mp.re(alpha)]+[mp.mpf(0)]*(L-1)+[mp.mpf(1)]
  opts[str(L)]=[{'mode':'original','Ablock':pair(original),'Bblock':pair(full)},{'mode':'all_lifts','Ablock':pair(full),'Bblock':pair(original)}]
  prepared.append({'block':block['block'],'original_group':block['original_group'],'lift_groups':lift_idx,'options':opts})
 # Fair choice: all originals for A, all lifts for B. Endpoints are set exact
 # in the future global GPU combination only after the balanced identity.
 H=np.asarray(c['H_coefficients_ascending'],float);fairA=[mp.mpf(1)];fairB=[mp.mpf(1)]
 for b in prepared:
  fairA=conv(fairA,[mp.mpf(x) for x in b['options'][str(L)][0]['Ablock']['float64']]);fairB=conv(fairB,[mp.mpf(x) for x in b['options'][str(L)][0]['Bblock']['float64']])
 fairA=np.asarray([float(x) for x in fairA]);fairB=np.asarray([float(x) for x in fairB]);fairA/=fairA[0];fairB/=fairB[0];HL=np.zeros((len(H)-1)*L+1);HL[::L]=H
 out={'format':'balanced-ds-block-options-v1','candidate':c['name'],'L':L,'mp_dps':100,'block_count':len(prepared),'blocks':prepared,'fair_cpu_highprec_H_error':float(abs(fairA-H).max()),'fair_cpu_highprec_HxL_error':float(abs(fairB-HL).max()),'fair_accept_lt_1e-9':bool(max(abs(fairA-H).max(),abs(fairB-HL).max())<1e-9),'mixed_radix_order':'block order; code least-significant digit is block0 option index','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'input_sha256':hashlib.sha256(a.input.read_bytes()).hexdigest(),'claim_boundary':'Prepared numerical coefficients only; GPU DS pass needs independent fair MPS calibration and CPU sample checks.'}
 a.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k]for k in ('candidate','L','mp_dps','block_count','fair_cpu_highprec_H_error','fair_cpu_highprec_HxL_error','fair_accept_lt_1e-9')},indent=2))
if __name__=='__main__':main()
