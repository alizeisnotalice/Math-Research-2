#!/usr/bin/env python3
"""Exact family: H binary and crossed factors positive, but C not binary."""
from __future__ import annotations
import hashlib,json
from pathlib import Path

HERE=Path(__file__).resolve().parent
def conv(a,b):return [sum(a[i]*b[k-i] for i in range(max(0,k-len(b)+1),min(k,len(a)-1)+1)) for k in range(len(a)+len(b)-1)]
def stretch(a,r):
 z=[0]*((len(a)-1)*r+1);z[::r]=a;return z
def poly(terms):
 z=[0]*(max(terms)+1)
 for i,c in terms.items():z[i]=c
 return z
def instance(p,q):
 F=poly({0:1,q-p:1,q:-1,q+p:1,2*q:1})
 G=poly({0:1,p:1,q:1,p+q:1})
 H=conv(F,G);A=conv(stretch(F,p),stretch(G,q));B=conv(stretch(G,p),stretch(F,q));C=conv(stretch(H,p),stretch(H,q))
 assert F[0]==G[0]==F[-1]==G[-1]==1 and F[q]==-1
 assert min(A)>=0 and min(B)>=0 and H[0]==H[-1]==1
 return {'p':p,'q':q,'F':F,'G':G,'H':H,'H_support':[i for i,c in enumerate(H) if c],
         'H_binary':all(c in (0,1) for c in H),'min_A':min(A),'min_B':min(B),
         'C_coefficient_at_pq':C[p*q],'C_binary':all(c in (0,1) for c in C)}
def main():
 target=instance(3,5)
 assert target['H_binary'] and target['C_coefficient_at_pq']==2 and not target['C_binary']
 assert target['H_support']==[0,2,3,5,7,8,10,11,13,15,16,18]
 primes=[3,5,7,11,13,17,19,23,29]
 scan=[]
 for i,p in enumerate(primes):
  for q in primes[i+1:]:
   z=instance(p,q)
   scan.append({'p':p,'q':q,'H_binary':z['H_binary'],'C_binary':z['C_binary']})
 out={'general_family':'F=1+x^(q-p)-x^q+x^(q+p)+x^(2q); G=(1+x^p)(1+x^q)',
      'target_p3_q5':target,'prime_scan_pq_le_29':scan,
      'exact_cross_positivity_explanation':{
       'A':'(1-x^(pq)+x^(2pq))*(1+x^(pq))=1+x^(3pq); remaining F terms are positive.',
       'B':'The middle q-fibre triplet x^(q(q-p))-x^(q^2)+x^(q(q+p)), after multiplication by 1+x^(pq), cancels both negative occurrences; endpoint terms remain nonnegative.'},
      'scope':'This is an exact counterexample to positivity descent after dropping C binary, while retaining H binary at (p,q)=(3,5). It is not a counterexample to P because C is not binary.'}
 path=HERE/'binary-H-nonbinary-C-sharpness.json';path.write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'target':{'H_binary':target['H_binary'],'C_binary':target['C_binary'],'C_pq':target['C_coefficient_at_pq'],'minA':target['min_A'],'minB':target['min_B']},'binary_H_pairs':[z for z in scan if z['H_binary']]},indent=2))
if __name__=='__main__':main()
