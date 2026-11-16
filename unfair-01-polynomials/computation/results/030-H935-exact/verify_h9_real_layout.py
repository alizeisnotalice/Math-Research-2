#!/usr/bin/env python3
"""No-SymPy independent integer audit of H9's 16 cycle masks and 80 layouts."""
from __future__ import annotations
import hashlib,json
from fractions import Fraction
from itertools import combinations
from math import gcd
from pathlib import Path

HERE=Path(__file__).resolve().parent; FULL=(1<<37)-1
FAMILIES=((16,2,18,3),(21,3,24,5),(29,3,32,5))

def trim(a):
 a=a[:]
 while len(a)>1 and a[-1]==0:a.pop()
 return a
def divexact(a,b):
 a=a[:];q=[0]*max(1,len(a)-len(b)+1)
 while len(a)>=len(b):
  t=a[-1]//b[-1];s=len(a)-len(b);q[s]=t
  for i,v in enumerate(b):a[s+i]-=t*v
  a=trim(a)
 assert a==[0]
 return trim(q)
def divisors(n):return [d for d in range(1,n+1) if n%d==0]
def cyclotomic(n,cache={}):
 if n in cache:return cache[n]
 a=[-1]+[0]*(n-1)+[1]
 for d in divisors(n):
  if d<n:a=divexact(a,cyclotomic(d))
 cache[n]=a;return a
def rem(a,m):
 a=trim(a)
 while len(a)>=len(m):
  t=a[-1];s=len(a)-len(m)
  for i,v in enumerate(m):a[s+i]-=t*v
  a=trim(a)
 return a+[0]*(len(m)-1-len(a))
def rref(a):
 a=[[Fraction(v) for v in row] for row in a];r=0;p=[];nr=len(a);nc=len(a[0])
 for c in range(nc):
  q=next((i for i in range(r,nr) if a[i][c]),None)
  if q is None:continue
  a[r],a[q]=a[q],a[r];f=a[r][c];a[r]=[v/f for v in a[r]]
  for i in range(nr):
   if i!=r:
    f=a[i][c]
    if f:a[i]=[v-f*w for v,w in zip(a[i],a[r])]
  p.append(c);r+=1
  if r==nr:break
 return a,p
def det(a):
 a=[row[:] for row in a];sign=1;prev=1
 for k in range(len(a)-1):
  if a[k][k]==0:
   q=next(i for i in range(k+1,len(a)) if a[i][k]);a[k],a[q]=a[q],a[k];sign=-sign
  pivot=a[k][k]
  for i in range(k+1,len(a)):
   for j in range(k+1,len(a)):a[i][j]=(a[i][j]*pivot-a[i][k]*a[k][j])//prev
  prev=pivot
  for i in range(k+1,len(a)):a[i][k]=0
 return sign*a[-1][-1]
def gword(eps,eta):
 z=0
 for e,t,(a,ac,b,bc) in zip(eps,eta,FAMILIES):
  if e:z|=((1<<ac)-1)<<a
  if t:z|=((1<<bc)-1)<<b
 return z
def decode(z):
 eps=[];eta=[]
 for a,ac,b,bc in FAMILIES:
  aa=(z>>a)&((1<<ac)-1);bb=(z>>b)&((1<<bc)-1)
  assert aa in (0,(1<<ac)-1) and bb in (0,(1<<bc)-1)
  eps.append(int(aa>0));eta.append(int(bb>0))
 return eps,eta

def main():
 groups=[]
 for n in (5,15,25):
  for j in range(1,(n-1)//2+1):
   if gcd(n,j)==1:groups.append((n,j,75*j//n))
 assert len(groups)==16
 phi75=cyclotomic(75);assert len(phi75)==41
 cols=[]
 for _,_,k in groups:
  v=[0]*76;v[k]+=1;v[75-k]+=1;cols.append(rem(v,phi75))
 M=[[cols[c][r] for c in range(16)] for r in range(40)]
 _,pivcols=rref(M);assert len(pivcols)==12
 sub=[[M[r][c] for c in pivcols] for r in range(40)]
 _,pivrows=rref(list(map(list,zip(*sub))));assert len(pivrows)==12
 minor=det([[M[r][c] for c in pivcols] for r in pivrows]);assert minor!=0
 blocks=[22,41,26176,39296]
 def matvec(mask):return [sum(M[r][c]*((mask>>c)&1) for c in range(16)) for r in range(40)]
 assert all(not any(matvec(z)) for z in blocks)
 supports=[[i for i in range(16) if z>>i&1] for z in blocks]
 assert sum(map(len,supports))==16 and not any(set(a)&set(b) for a,b in combinations(supports,2))
 zero=sorted(sum(blocks[i] for i in range(4) if q>>i&1) for q in range(16))
 assert all(not any(matvec(z)) for z in zero)
 index={k:i for i,(_,_,k) in enumerate(groups)};want={(1,3):22,(2,3):41,(1,5):26176,(2,5):39296};binding={}
 for (j,p),w in want.items():
  roots=[q for q in range(1,75) if (p*q-15*j)%75==0]
  got=sum(1<<index[min(q,75-q)] for q in roots);assert got==w;binding[f'j{j}_p{p}']=got
 raw=json.loads((HERE/'H9-real-necessary-layouts.json').read_text())['canonical_layouts'];actual={r['mask_u37'] for r in raw};assert len(actual)==80
 expected=set(zero)
 for c in zero:
  for bits in range(4):
   eps=[1,bits&1,(bits>>1)&1];eta=[0,1-eps[1],1-eps[2]];expected.add(c|gword(eps,eta))
 assert actual==expected
 for r in raw:
  z=r['mask_u37'];eps,eta=decode(z);assert z&0xffff in zero and r['selected_group_indices']==[i for i in range(37) if z>>i&1] and r['complement_mask_u37']==z^FULL
  if r['class_kind']=='e0_e2':assert eps==eta==[0,0,0]
  else:assert eps[0]==1 and eta[0]==0 and all(a+b==1 for a,b in zip(eps,eta))
 standard=63|gword([1,1,1],[0,0,0]);assert standard in actual and 0 in actual
 out={'standard_library_only':True,'phi75_degree':40,'matrix_rank':12,'minor_rows':pivrows,'minor_columns':pivcols,'minor_determinant':minor,'kernel_block_masks':blocks,'Phi5_lift_bindings':binding,'zero_cycle_boolean_masks':zero,'exact_expected_80_set_equals_ledger':True,'constant_pair_present':True,'standard_pair_present':True,'claim_boundary':'Necessary-layout audit only; no coefficient sign test.'}
 p=HERE/'H9-real-layout-independent-verification.json';p.write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({**out,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()},indent=2))
if __name__=='__main__':main()
