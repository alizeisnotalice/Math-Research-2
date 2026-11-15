"""Exact finite checks supporting the analytic stretch lemmas (not a search)."""
from itertools import permutations
from math import comb
import json
import sympy as S

def check_group(G,d):
 for i in range(d):
  stab=[g for g in G if g[i]==i];h=len(stab);assert h%(d-1)==0;c=h//(d-1)
  for j in range(d):
   for k in range(d):
    actual=sum(g[j]==k for g in stab)
    expected=h if j==k==i else (c if j!=i and k!=i else 0)
    assert actual==expected
 return {'degree':d,'order':len(G),'all_stabilizer_matrix_entries_verified':True}
A4={p for p in permutations(range(4)) if sum(p[i]>p[j] for i in range(4) for j in range(i+1,4))%2==0}
AGL={tuple((a*i+b)%5 for i in range(5)) for a in range(1,5) for b in range(5)}
counts=[]
for L in (1,3,5):
 total=0;blockbits=L+1;liftmask=(1<<L)-1
 for mask in range(1<<(3*blockbits)):
  es=[]
  for i in range(3):
   block=(mask>>(i*blockbits))&((1<<blockbits)-1)
   es.append(L*(block&1)+((block>>1)&liftmask).bit_count())
  if len(set(es))==1:total+=1
 formula=2*sum(comb(L,e)**3 for e in range(1,L))+10
 assert total==formula
 counts.append({'L':L,'all_raw_masks':1<<(3*blockbits),'balanced_ordered':total,'nonconstant_complement_representatives':(total-2)//2})
y=S.symbols('y');Hminus=y**6-y**5+y**2-y+1
assert S.expand(Hminus-(2*y*y-2*y+1)-y*(1-y)*(1-y**4))==0
assert S.expand(Hminus-y**5*(y-1)-(y*y-y+1))==0
large=[{'L':L,'nonconstant_complement_representatives':sum(comb(L,e)**3 for e in range(1,L))+4} for L in (7,8,9,10,11,12)]
print(json.dumps({'stabilizer_checks':[check_group(A4,4),check_group(AGL,5)],'small_complete_integer_mask_checks':counts,'large_exact_counts':large,'H_no_real_roots_decomposition_verified':True,'scope':'Group matrices, identities and combinatorial counts; no GPU non-hit is an exclusion'},indent=2))
