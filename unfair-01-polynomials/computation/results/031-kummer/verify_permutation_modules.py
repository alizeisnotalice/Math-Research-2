#!/usr/bin/env python3
"""Small exact F_p permutation-module checks used in the Kummer note."""
from itertools import product, permutations
import hashlib,json
from pathlib import Path

HERE=Path(__file__).resolve().parent
def rank(rows,p):
 a=[list(r) for r in rows if any(x%p for x in r)];r=c=0
 while r<len(a) and c<len(a[0]):
  q=next((i for i in range(r,len(a)) if a[i][c]%p),None)
  if q is None:c+=1;continue
  a[r],a[q]=a[q],a[r];inv=pow(a[r][c]%p,-1,p);a[r]=[(z*inv)%p for z in a[r]]
  for i in range(len(a)):
   if i!=r and a[i][c]%p:
    f=a[i][c]%p;a[i]=[(u-f*v)%p for u,v in zip(a[i],a[r])]
  r+=1;c+=1
 return r
def orbit_span(v,p):
 return rank([tuple(v[i] for i in perm) for perm in permutations(range(len(v)))],p)
def main():
 # p not dividing n: every nonzero sum-zero vector is nonconstant and has U orbit.
 n,p=5,3;u=[v for v in product(range(p),repeat=n) if sum(v)%p==0]
 assert all(orbit_span(v,p)==n-1 for v in u if any(v))
 # p dividing n: constants are inside U and form the exceptional T submodule.
 n2,p2=5,5;u2=[v for v in product(range(p2),repeat=n2) if sum(v)%p2==0]
 constants={tuple([a]*n2) for a in range(p2)}
 assert all(orbit_span(v,p2)==1 for v in constants if any(v))
 assert all(orbit_span(v,p2)==n2-1 for v in u2 if v not in constants)
 out={'p_not_dividing_n_example':{'n':n,'p':p,'U_size':len(u),'every_nonzero_U_vector_orbit_rank':n-1},
      'p_dividing_n_example':{'n':n2,'p':p2,'U_size':len(u2),'constant_T_size':len(constants),'nonzero_T_orbit_rank':1,'every_nonconstant_U_vector_orbit_rank':n2-1},
      'scope':'Finite exact illustrations of the general transposition-difference proof, not a proof of Kummer non-pth-power hypotheses.'}
 path=HERE/'permutation-module-exact-check.json';path.write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),**out},indent=2))
if __name__=='__main__':main()
