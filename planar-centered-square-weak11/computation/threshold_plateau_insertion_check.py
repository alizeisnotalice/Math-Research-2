#!/usr/bin/env python3
"""Deterministic exact check of the single-scale threshold-plateau example."""
from fractions import Fraction as F
from pathlib import Path
import argparse, json, platform, time

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--out',required=True); args=ap.parse_args()
    start=time.perf_counter(); epsilon=F(1,1024); r=F(1,2)
    support_f=(-F(1),F(1)); support_g=(-F(1,4),F(1,4))
    assert 0<epsilon<=F(1,4)
    assert 4*epsilon<=1
    area_f=(support_f[1]-support_f[0])**2
    area_g=(support_g[1]-support_g[0])**2
    mass_g=4*epsilon*area_g
    qarea=4*r*r
    assert mass_g==epsilon and qarea==1
    corner_checks=[]
    for u in support_g:
        for v in support_g:
            for a in (u,v):
                assert support_f[0]<=a-r<=support_g[0]
                assert support_g[1]<=a+r<=support_f[1]
            corner_checks.append([str(u),str(v)])
    average=(qarea+mass_g)/qarea
    assert average==1+epsilon>1
    result=dict(status='passed',python=platform.python_version(),method='Fraction exact arithmetic; deterministic single input',
      epsilon=str(epsilon),f_square=[str(a) for a in support_f],g_square=[str(a) for a in support_g],
      f_density='1',g_density=str(4*epsilon),f_mass=str(area_f),g_mass=str(mass_g),
      radius=str(r),square_area=str(qarea),receiver_area_lower=str(area_g),
      mixed_average=str(average),new_area_per_mass_lower=str(area_g/epsilon),
      checked_receiver_corners=corner_checks,random_seed=None,tolerance=None,
      elapsed_seconds=time.perf_counter()-start,
      scope='Exact scalar and endpoint checks; the all-epsilon and all-center argument is analytic in threshold-plateau-insertion.md')
    Path(args.out).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__': main()
