#!/usr/bin/env python3
"""Exact maximum of a symmetric rational quadratic form on a simplex.
Enumerate all face stationary points with nonsingular bordered KKT systems.
A minimal-support global maximizer has such a system; see PROOF.md.
"""
from fractions import Fraction as F
from pathlib import Path
from itertools import combinations
import json,sys,hashlib

def solve(A,b):
 a=[list(map(F,r))+[F(y)] for r,y in zip(A,b)];n=len(a)
 for i in range(n):
  p=next((p for p in range(i,n) if a[p][i]),None)
  if p is None:return None
  a[p],a[i]=a[i],a[p];z=a[i][i];a[i]=[t/z for t in a[i]]
  for j in range(n):
   if j!=i:
    z=a[j][i];a[j]=[x-z*y for x,y in zip(a[j],a[i])]
 return [row[-1] for row in a]
def audit(path):
 x=json.loads(Path(path).read_text());n=x['n'];N=n*n;A=int(x['A2'])
 if 'covariance_numerators' in x:
  B=[[F(int(z),int(x['covariance_common_denominator'])) for z in row] for row in x['covariance_numerators']]
 else:
  l=list(map(int,x['cut_marginal_numerators']));h=[[int(z) for z in row] for row in x['ordered_pair_cross_numerators']]
  B=[[F(2*A*h[j][k]-l[j]*l[k],4*A*A) for k in range(n-1)] for j in range(n-1)]
 assert len(B)==n-1 and B==[list(row)for row in zip(*B)]
 candidates=[];singular=[];infeasible=[]
 for m in range(1,n):
  for inds in combinations(range(n-1),m):
   K=[[B[j][k]for k in inds]+[F(-1)] for j in inds]+[[F(1)]*m+[F(0)]]
   sol=solve(K,[0]*m+[1])
   if sol is None:singular.append(list(inds));continue
   w=sol[:-1];lam=sol[-1]
   if min(w)<0:infeasible.append(list(inds));continue
   assert sum(w)==1
   val=sum(w[j]*w[k]*B[inds[j]][inds[k]] for j in range(m) for k in range(m));assert val==lam
   candidates.append({'face':list(inds),'weights':[str(z)for z in w],'value':str(val)})
 best=max(candidates,key=lambda x:F(x['value']))
 return {'n':n,'input_sha256':hashlib.sha256(Path(path).read_bytes()).hexdigest(),'faces':2**(n-1)-1,'singular_faces':singular,'infeasible_stationary_faces':infeasible,'candidates':candidates,'maximum':best,'negative_copositive':F(best['value'])<0,'maximum_decimal':float(F(best['value']))}
if __name__=='__main__':
 out={'arithmetic':'fractions.Fraction exact','results':[audit(p)for p in sys.argv[1:]]};print(json.dumps(out,indent=2))
