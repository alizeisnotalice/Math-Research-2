#!/usr/bin/env python3
"""Exact singleton tests for linear single-column potential obstruction."""
from fractions import Fraction as F
from pathlib import Path
import json,platform,hashlib

def width(a,b,c,d):return max(F(0),min(b,d)-max(a,c))
def main():
  rows=[];tests=0
  for n in range(2,10):
    h=F(1,2**n);T=F(2**n);R=T;step=F(1,2**(4*n));N=2**(n-2)
    H=F(0);seen=set()
    for l in range(1,N+1):
      r=l*h;m=4*r*r
      assert R>r and T>2*r and m<=F(1,4)
      assert (m/step).denominator==1
      for j in [l,1-l]:
        assert j not in seen;seen.add(j)
        # Two distinct periodic lifts at y=0,T supply the lower bound.
        area=sum((width(j*h-r,j*h+r,F(0),h)*width(y-r,y+r,F(0),T) for y in [F(0),T]),F(0))
        assert area==2*l*h*h
        assert area/m==F(1,2*l)
        # A strict interior midpoint has a witnessing radius strictly below r.
        x=h/2;y=r/2
        distance=max(abs(x-j*h),abs(y))
        assert distance<r
        rw=(distance+r)/2
        assert 0<rw<r<R and m>4*rw*rw
        tests+=1
      H+=F(1,l)
    assert len(seen)==2*N
    rows.append({'n':n,'N':N,'harmonic_lower_bound_exact':str(H),'maximum_test_mass':str(4*(N*h)**2),'tested_indices_min':min(seen),'tested_indices_max':max(seen)})
  k=4;z=F(5,3)
  alpha={j:F(j*j-3*j,7) for j in range(-k,k)}
  beta={j:alpha.get(j-1,F(0))-alpha.get(j,F(0))+(z if j==0 else 0) for j in range(-k,k+1)}
  assert sum(beta.values(),F(0))==z
  assert beta[-k]==-alpha[-k] and beta[k]==alpha[k-1]
  out={'python':platform.python_version(),'single_process':True,'arithmetic':'fractions.Fraction','random_seed':None,'tolerance':None,'singleton_tests':tests,'rows':rows,'telescoping_beta_sum':str(z),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
  dest=Path(__file__).parent/'results/singleton_column_potential_obstruction_exact.json';dest.parent.mkdir(parents=True,exist_ok=True);dest.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
