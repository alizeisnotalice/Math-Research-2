#!/usr/bin/env python3
"""Single-core exact spot checks for fixed-physical-range obstruction."""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import json,hashlib,platform

def length(a,b,c,d):return max(F(0),min(b,d)-max(a,c))
def main():
 rows=[];subsets_checked=0;averages_checked=0
 for L in [1,3]:
  D=F(2*L-1,2);assert L>D
  for p in range(3):
   m=16*L*L*4**p;n=(m-1).bit_length();h=F(1,2**n);T=F(2**n);R=T;q=L*2**n
   assert F(m)/F(1,2**(4*n))==m*2**(4*n)
   rmax=F(L*4**p);assert 4*rmax<m<=T and 2*rmax<R
   lower=F(0)
   for s in range(p+1):
    r=F(L*2**(p+s));M=2**(p+s);K=4**s;t=2*r
    assert K<=M and K*m==4*t*t
    assert q*M==r/h
    averaged=(F(M,K)*4*r*h/m);assert averaged==F(1,4*q)
    lower+=q*averaged
    for residue in [0,q//2,q-1]:
     xs=[r+residue*h+j*L for j in range(M)]
     assert all(r<=x<2*r for x in xs)
     total=F(0);number=0
     weights=[F((i+1)**2,17) for i in range(M)]
     for indices in combinations(range(M),K):
      selected=[xs[i] for i in indices]
      assert all(b-a>=L>D for a,b in zip(selected,selected[1:]))
      assert all(0<x<R for x in selected)
      xa=max(selected)-t;xb=min(selected)+t
      area=sum((length(xa,xb,F(0),h)*length(y-t,y+t,F(0),T) for y in [F(0),T]),F(0))
      assert area==4*r*h
      x=h/2;y=r
      distance=max([abs(v-x) for v in selected]+[abs(y)])
      rw=(distance+t)/2
      assert distance<rw<t<R and K*m>4*rw*rw
      total+=sum((weights[i] for i in indices),F(0));number+=1;subsets_checked+=1
     assert total/number==F(K,M)*sum(weights,F(0));averages_checked+=1
   assert lower==F(p+1,4)
   rows.append({'L':L,'p':p,'n':n,'mass':m,'q':q,'annular_beta_lower_bound':str(lower)})
 # Explicit all-n choice n0=4+2ceil(log2 L).
 for L in [1,2,3,5,16]:
  n0=4+2*(L-1).bit_length()
  for n in range(n0,n0+9):
   p=(n-n0)//2
   assert 16*L*L*4**p<=2**n
 # Normalize arbitrary slot constants before forming alpha at fixed test mass.
 k=3;mass=F(7);z=F(5,2)
 phi0={j:F(j*j+1,3) for j in range(-k,k)}
 alpha={j:F(j+4,5)*mass for j in range(-k,k)}
 beta={j:alpha.get(j-1,F(0))-alpha.get(j,F(0))+(z if j==0 else 0) for j in range(-k,k+1)}
 assert sum(beta.values(),F(0))==z
 for position in range(-k,k+1):
  def potential(shift):
   return sum((phi0[j]+F(j+4,5)*(mass if j+shift==position else F(0))**2 for j in range(-k,k)),F(0))
  delta=potential(1)-potential(0)
  assert delta==mass*(beta[position]-(z if position==0 else 0))
 out={'python':platform.python_version(),'single_process':True,'arithmetic':'fractions.Fraction','random_seed':None,'tolerance':None,'subset_geometry_checks':subsets_checked,'combination_average_checks':averages_checked,'layers':rows,'all_n_mass_choices_checked':45,'slot_constant_normalization_checks':2*k+1,'telescoping_beta_sum':str(z),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 dest=Path(__file__).parent/'results/fixed_physical_range_potential_obstruction_exact.json';dest.parent.mkdir(parents=True,exist_ok=True);dest.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
