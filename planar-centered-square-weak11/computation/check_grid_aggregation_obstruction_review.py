#!/usr/bin/env python3
"""Independent exact checks of fixed-support grid aggregation obstruction."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,platform

def main():
    rows=[];pairs=0
    assert F(10,7)**2>2 # exact certificate sqrt(2)<10/7
    gap=F(1,14)-F(3,1024)
    assert gap>F(1,16)
    for n in range(3,9):
      h=F(1,2**n);s=h/8;eps=s/4096
      K=F(32)/(h*h);assert K.denominator==1
      mass=2*s*s;assert K*mass==1
      rho=mass/(h*h);assert rho==F(1,32)
      # Every neighboring cell pair, all four source choices; farther cells are separated by >=h.
      for dx in [-1,0,1]:
       for dy in [-1,0,1]:
        if dx==dy==0: continue
        for a in [-1,1]:
         for b in [-1,1]:
          dist=max(abs(dx*h+(b-a)*s/2),abs(dy*h))
          assert dist>=h-s
          pairs+=1
      r0=(h-s)/2;delta=s/2
      tail=rho*(1+(h+2*delta)/(2*r0))**2
      assert tail==F(8,49)<1
      aggregated_tail=rho*(1+h/(2*(h/2)))**2
      assert aggregated_tail==F(1,8)<1
      step=F(1,2**(4*n));units=mass/step
      assert units.denominator==1 and units==2**(2*n-5)
      assert s/2+eps<h/2 and eps<s/4
      # Exact rational enclosure sqrt2<10/7 certifies pair rectangle loss <=8s eps.
      assert 4*(2*F(10,7)-1)<8
      loss=24*K*s*eps;assert loss==F(3,1024)
      # Singletons after aggregation have radius sqrt(mass)/2 < h/2.
      assert mass<h*h
      # Uniform density observation, not required by the original theorem.
      density=s*s/(4*eps*eps);assert density==2**22
      rows.append({'n':n,'h':str(h),'K':int(K),'cell_mass':str(mass),'quantization_units':int(units),'cross_cell_tail':str(tail),'aggregated_tail':str(aggregated_tail),'smoothing_loss_bound':str(loss)})
    out={'python':platform.python_version(),'single_process':True,'arithmetic':'fractions.Fraction','random_seed':None,'tolerance':None,'neighbor_source_pairs_checked':pairs,'layers':rows,'strict_gap_lower_bound':str(gap),'excess_over_1over16':str(gap-F(1,16)),'uniform_Linfty':2**22,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    dest=Path(__file__).parent/'results/grid_aggregation_obstruction_review_exact.json';dest.parent.mkdir(exist_ok=True,parents=True);dest.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
