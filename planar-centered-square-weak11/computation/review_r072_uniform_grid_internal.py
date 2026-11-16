#!/usr/bin/env python3
"""Fraction shape/cap controls for the uniform-grid analytic proof."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from math import isqrt
from pathlib import Path
import sys
import time


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();start=time.perf_counter();spec=json.loads(a.input.read_text());cases=[];total=0
    den=1<<spec['sqrt_bits'];q=isqrt(2*den*den);rootlo,roothi=F(q,den),F(q+1,den);assert rootlo**2<=2<roothi**2
    for n in spec['grid_sizes']:
        assert isinstance(n,int) and n>=2
        w=1-F(1,n*n);h=F(n+1,n*n);t=F(n,n+1)
        assert w/(n-1)==h
        eligible=adjacent=boundtests=fcaptests=0
        for k in range(1,n+1):
            ak=(k-1)*h+F(1,n*n)
            assert 0<ak<=1
            for l in range(1,n+1):
                total+=1
                mass=F(k*l,n*n);diam=max(k-1,l-1)*h
                direct=mass>diam**2
                expected=k==l or (abs(k-l)==1 and min(k,l)*(2*n+1)<n*n)
                assert direct==expected
                if direct:
                    eligible+=1;adjacent+=k!=l
                    # sqrt(k*l)<=A+sqrt(2), isolating the only radical.
                    A=2*n*(k-1)*h;z=k*l-A*A-2
                    assert A>=0 and (z<=0 or z*z<=8*A*A)
                    boundtests+=1
                al=(l-1)*h+F(1,n*n);f=mass-t*ak*al
                assert f<=1-t
                if k<n:
                    nextf=F((k+1)*l,n*n)-t*(ak+h)*al
                    assert nextf-f==(1-F(l,n))/(n*n)
                fcaptests+=1
        assert (n-1)*h+F(1,n*n)==1
        assert F(n*n,n*n)-t==1-t
        if n>=3:assert F(2,n*n)>h*h # exact maximizer of one-face extension
        support=(1+w)**2
        cases.append(dict(n=n,all_ordered_shapes=n*n,eligible_shapes=eligible,adjacent_shapes=adjacent,extension_radical_comparisons=boundtests,diffuse_scalar_cap_pairs=fcaptests,internal_halo_containing_square_area_lower=str((w+rootlo/n)**2),internal_halo_containing_square_area_upper=str((w+roothi/n)**2),internal_halo_containing_square_upper_display_only=float((w+roothi/n)**2),support_carrier_area=str(support),support_carrier_minus_7_over_2=str(support-F(7,2)),exact_maximum_one_face_extension='sqrt(2)/(2*n)' if n>=3 else '1/(2*n)'))
    result=dict(status='PASS_EXACT_UNIFORM_GRID_SHAPE_AND_LOCALIZATION_CONTROLS',source_sha256=sha(Path(__file__)),input_sha256=sha(a.input),python=sys.version,ordered_shapes_checked=total,cases=cases,elapsed_seconds=time.perf_counter()-start,scope='Exact k,l eligibility and directed scalar localization checks; general statements require accompanying analytic proof. Diffuse-family scalar cap identities are not a finite-atomic residual-halo transfer, no actual general Delta or C2 bound.')
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(status=result['status'],ordered_shapes_checked=total,seconds=result['elapsed_seconds'])))


if __name__=='__main__':main()
