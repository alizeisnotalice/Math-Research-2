#!/usr/bin/env python3
"""Exact F2^4 plateaued switch vector with signed cancellation on a true difference."""
import json, platform, sys, hashlib, time
from pathlib import Path
N,q=16,4
def fwht(a):
 a=list(a);d=1
 while d<N:
  for i in range(0,N,2*d):
   for j in range(i,i+d):a[j],a[j+d]=a[j]+a[j+d],a[j]-a[j+d]
  d*=2
 return a
def main():
 started=time.monotonic()
 # A={x_3=0}; h=1_A (-1)^(x_1 x_2).
 A={x for x in range(N) if not ((x>>3)&1)}
 h=[0]*N
 for x in A: h[x]= -1 if (((x>>1)&1) and ((x>>2)&1)) else 1
 f=[(((x>>1)&1)*((x>>2)&1)) ^ ((x&1)*((x>>3)&1)) for x in range(N)]
 g=[f[x] ^ int(x in A) for x in range(N)]
 assert all(h[x]==(1-2*f[x] if x in A else 0) for x in range(N))
 assert set(map(abs,fwht([1-2*z for z in f])))=={q}
 assert set(map(abs,fwht([1-2*z for z in g])))=={q}
 W=fwht(h); S={z for z,w in enumerate(W) if w}
 assert len(A)==len(S)==8 and {abs(w) for w in W if w}=={4}
 delta=2 # e_1, a genuine difference inside A+A
 signed=sum(h[x]*h[x^delta] for x in range(N))
 unsigned=sum((x in A) and ((x^delta) in A) for x in range(N))
 diff={x^y for x in A for y in A}
 assert signed==0 and unsigned==8 and delta in diff
 # Kernel <e_3> avoids A+A, so the example is not a counterexample to projectability.
 assert 8 not in diff
 out={'base_and_flip_bent_verified':True,'base_anf':'x1*x2+x0*x3','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-started,'python':sys.version,'platform':platform.platform(),'support_size':len(A),'walsh_support_size':len(S),'nonzero_walsh_abs':sorted({abs(w) for w in W if w}),'delta':delta,'delta_in_A_plus_A':delta in diff,'signed_autocorrelation':signed,'unsigned_pair_count':unsigned,'projecting_kernel_generator':8,'kernel_avoids_A_plus_A':8 not in diff}
 json.dump(out,Path(__file__).with_name('signed_autocorrelation_obstacle_k2.json').open('w'),indent=2,sort_keys=True);print(json.dumps(out,indent=2,sort_keys=True))
if __name__=='__main__':main()
