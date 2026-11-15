#!/usr/bin/env python3
"""Exact arithmetic support for R53 radius-split contraction."""
from fractions import Fraction as F
from pathlib import Path
import argparse,hashlib,json,platform,sys,time

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();start=time.perf_counter()
 B=F(1994727599733432151936139,604462909807314587353088)
 delta=F(1,1250);r0=F(1,25);n=9;factor=(1+delta/r0)**2
 upper=F(33,10)*factor+4*n*(r0+delta)**2
 assert factor==F(2601,2500) and upper==F(10916397,3125000)
 assert F(7,2)-upper==F(21103,3125000) and B*factor+4*n*(r0+delta)**2<upper<F(7,2)
 checks=0
 for e in [F(0),F(1,1250),F(1,7),F(2)]:
  for s in [F(1,25),F(1,3),F(1),F(5)]:
   for r in [s,2*s,10*s,100*s]:
    assert (1+e/s)**2*r*r >= (r+e)**2;checks+=1
   # Derivative of (s+e)^2*(B/s^2+4n), algebraically factored.
   lhs=2*(s+e)*(B/s**2+4*n)-2*B*(s+e)**2/s**3
   rhs=2*(s+e)*(4*n-B*e/s**3)
   assert lhs==rhs;checks+=1
 result=dict(status='PASS',python=sys.version,platform=platform.platform(),arithmetic='fractions.Fraction exact; deterministic; no tolerances',checks=checks,delta_ratio=str(delta),r0_ratio=str(r0),block_count=n,mass_inflation=str(factor),upper_exact=str(upper),margin_below_7_2=str(F(7,2)-upper),strict_bound_using_certificate=str(B*factor+4*n*(r0+delta)**2),elapsed_seconds=time.perf_counter()-start,source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__).resolve(),Path(__file__).resolve().parents[1]/'proofs/r053-radius-split-contraction.md']},scope='Scalar verification supporting analytic radius-split theorem and conditional nine-tight-block corollary, not general C2.')
 a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
