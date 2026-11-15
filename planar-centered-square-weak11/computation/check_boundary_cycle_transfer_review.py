#!/usr/bin/env python3
"""Independent exact checks for transfer/rounding/potential; not general proof."""
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import json,platform,hashlib

def main():
    capchecks=0
    for m in [F(1,16),F(1),F(2),F(4),F(5),F(9)]:
      for R in [F(1,2),F(1),F(3,2)]:
       for D in range(1,65):
        s2=m/4
        v=F(D*D)*s2
        fl=isqrt(v.numerator//v.denominator)
        t=min(R,max(F(0),F(fl-1,D)))
        assert t>=0 and t<=R and t*t<=s2
        if t>0:
          assert (t+F(1,D))**2<=s2
          assert (t+F(1,2*D))**2<s2
          capchecks+=1
    rho=F(13,10);r=F(1,8)
    points=[(x,y) for x in [F(1,4),F(1,2),F(3,4)] for y in [F(1,4),F(1,2),F(3,4)]]
    quant=[]
    for n in range(3,8):
      h=F(1,2**n);step=F(1,2**(4*n));count=2**(2*n)
      oldcell=rho*h*h;newcell=(oldcell//step)*step
      mass=count*newcell;loss=rho-mass
      assert F(0)<=loss<count*step
      assert mass>0 and newcell<=2**n
      rad=r+h;old_witness_mass=rho*4*r*r
      slack=old_witness_mass-4*r*r
      sufficient=loss+8*r*h+4*h*h<slack
      avgs=[]
      for x,y in points:
       # Count source grid points in the enlarged closed witness directly.
       nx=sum(1 for j in range(2**n) if abs(j*h-x)<=rad)
       ny=sum(1 for j in range(2**n) if abs(j*h-y)<=rad)
       qmass=nx*ny*newcell
       if sufficient: assert qmass>4*rad*rad
       avgs.append(qmass/(4*rad*rad))
      assert all(a>1 for a in avgs)
      quant.append({'n':n,'mass':str(mass),'loss':str(loss),'uniform_margin_bound_sufficient':sufficient,'minimum_actual_average':str(min(avgs))})
    # Finite graph with a zero-cost positive-reward EDGE and zero cycle.
    edges=[(0,1,F(2),F(0)),(1,0,F(0),F(2)),(1,1,F(1,2),F(1)),(1,2,F(0),F(0)),(2,1,F(0),F(1)),(0,0,F(0),F(0))]
    z=F(1);pot=[F(0),F(0),F(0)]
    def dfs(u,seen,value):
      pot[u]=max(pot[u],value)
      for a,b,reward,cost in edges:
       if a==u and b not in seen: dfs(b,seen|{b},value+reward-z*cost)
    for start in range(3):dfs(start,{start},F(0))
    assert all(a-z*b<=pot[v]-pot[u] for u,v,a,b in edges)
    assert F(2)-F(3,4)*2>0 # positive cycle forbids potential at z=3/4
    out={'python':platform.python_version(),'single_process':True,'arithmetic':'Fraction/isqrt','random_seed':None,'tolerance':None,'positive_inner_radius_checks':capchecks,'quantization':quant,'graph_potential_at_z_1':[str(v) for v in pot],'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    dest=Path(__file__).parent/'results/boundary_cycle_transfer_review_exact.json';dest.parent.mkdir(exist_ok=True,parents=True);dest.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
