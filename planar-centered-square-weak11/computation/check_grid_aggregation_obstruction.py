#!/usr/bin/env python3
"""Exact finite checks for the uniform grid aggregation obstruction.
No global halo enumeration is performed; the all-n proof is analytic.
"""
from fractions import Fraction as F
from pathlib import Path
from itertools import product
import hashlib,json,os,platform,sys

def main():
    assert os.path.ismount('/Volumes/data')
    checks=[]; pairs=0
    for n in range(3,9):
        h=F(1,2**n); s=h/8; eps=s/4096
        K=(8/h)*(4/h); assert K.denominator==1
        mass=2*s*s; delta=s/2; r0=(h-s)/2; rho=mass/(h*h)
        assert K*mass==1 and rho==F(1,32)
        bound=rho*(1+(h+2*delta)/(2*r0))**2
        assert bound==F(8,49)<1
        quant_tail=rho*(1+h/(2*(h/2)))**2
        assert quant_tail==F(1,8)<1
        # A center and all eight neighboring cells; farther cells have larger
        # separation by an explicit coordinate bound in the analytic proof.
        for di,dj in product((-1,0,1),repeat=2):
            if not(di or dj): continue
            for a,b in product((-1,1),repeat=2):
                dist=max(abs(di*h+(b-a)*s/2),abs(dj*h))
                assert dist>=h-s
                pairs+=1
        d=F(1,2**(4*n)); multiple=mass/d
        assert multiple.denominator==1 and multiple==2**(2*n-5)
        assert s/2+eps<h/2 and eps<h/2
        total_loss=24*K*s*eps
        assert total_loss==F(3,1024)
        checks.append({'n':n,'h':str(h),'s':str(s),'epsilon':str(eps),'cells':int(K),
          'mass_per_cell':str(mass),'total_mass':str(K*mass),'quantization_multiple':int(multiple),
          'cross_cell_tail':str(bound),'quantized_tail':str(quant_tail),'area_loss_bound':str(total_loss)})
    # sqrt(2)<10/7, so beta-1-loss>1/14-3/1024>1/16.
    assert F(10,7)**2>2
    gap_lower=F(1,14)-F(3,1024)
    assert gap_lower>F(1,16)
    source=Path(__file__).resolve()
    result={'status':'PASS','python_version':sys.version,'platform':platform.platform(),
      'script_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
      'arithmetic':'fractions.Fraction','seed':None,'tolerance':None,
      'levels':checks,'neighbor_source_pairs_checked':pairs,
      'uniform_halo_loss_strict_lower_bound':'1/16','rational_intermediate_gap_lower_bound':str(gap_lower),
      'scope':'finite scalar and neighbor checks; all-n full halo and L1 claims use the analytic proof',
      'not_claimed':'no lower bound on C2-L_n; no large source-box enumeration'}
    out=source.parent/'results/recursive-boundary/grid-aggregation-obstruction.json'
    out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'levels':len(checks),'neighbor_source_pairs_checked':pairs,
      'cross_cell_tail':'8/49','quantized_tail':'1/8','uniform_halo_loss_strict_lower_bound':'1/16'},indent=2))
    print(out)
if __name__=='__main__': main()
