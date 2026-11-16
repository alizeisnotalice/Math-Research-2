#!/usr/bin/env python3
"""Exact all-center unit-cap certificate for the scaled R75 embedding.

No full halo union is recomputed. Root/compute own that distinct certificate.
"""
from fractions import Fraction as F
from itertools import product
import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys
import time


def axis_atoms(edges):
    edges=sorted(set(edges));out=[]
    for i,x in enumerate(edges):
        out.append((x,x,x))
        if i+1<len(edges):out.append((x,edges[i+1],(x+edges[i+1])/2))
    return out


def contains(rect,p):return rect[0]<=p[0]<=rect[1]and rect[2]<=p[1]<=rect[3]


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--input',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();start=time.perf_counter();raw=args.input.read_bytes();d=json.loads(raw)
    n=d['n'];assert n==3;e=F(1,n*n);h=F(n+1,n*n)
    pts=[(i*h,j*h)for i in range(n)for j in range(n)]+[tuple(map(F,p))for p in d['outside_positions']]
    masses=[e]*(n*n)+list(map(F,d['outside_masses']))
    axes=[axis_atoms([p[a]+s*F(1,2)for p in pts for s in (-1,1)])for a in (0,1)]
    cells=[];maximum=F(0)
    for ix,xa in enumerate(axes[0]):
        for iy,ya in enumerate(axes[1]):
            x,y=xa[2],ya[2];Q=(x-F(1,2),x+F(1,2),y-F(1,2),y+F(1,2))
            ids=[i for i,p in enumerate(pts)if contains(Q,p)]
            m=sum((masses[i]for i in ids),F(0));assert m<=1;maximum=max(maximum,m)
            cells.append({'axis_atom_indices':[ix,iy],'representative':[str(x),str(y)],'captured_indices':ids,'mass':str(m),'cap_slack':str(1-m)})
    assert len(cells)==225 and maximum==1
    Q=tuple(map(F,d['nonunit_test_square']));assert Q[1]-Q[0]==Q[3]-Q[2]==F(10,9)
    ids=[i for i,p in enumerate(pts)if contains(Q,p)]
    m=sum((masses[i]for i in ids),F(0));cost=(Q[1]-Q[0])**2
    assert m==F(4,3)and cost==F(100,81)and m-cost==F(8,81)>0
    l2,l5,u2,u5=map(F,[d['sqrt2_lower'],d['sqrt5_lower'],d['sqrt2_upper'],d['sqrt5_upper']])
    assert l2*l2<2<u2*u2 and l5*l5<5<u5*u5
    whole_lower=4*l2+l5/2-F(25,4)
    assert whole_lower==F(311,612)and whole_lower-F(1,2)==F(5,612)>0
    # Analytic clipped THREE-atom formula, not a full-S3 halo calculation.
    clipped_upper=4*u2/9+u5/54-F(209,324)
    assert clipped_upper<F(1,18)
    out={'schema':'r075-scaled-unit-qualification-result-v1','all_passed':True,
         'positive_source_count':len(pts),'positions':[[str(v)for v in p]for p in pts],'masses':list(map(str,masses)),
         'unit_center_axis_atoms':[[[str(v)for v in a]for a in aa]for aa in axes],
         'all_unit_center_cells':cells,'cell_count':len(cells),'maximum_unit_square_mass':str(maximum),
         'outside_center_hull':'one coordinate too far from every source, so no source captured; mass zero',
         'nonunit_violation':{'closed_square':list(map(str,Q)),'captured_indices':ids,'mass':str(m),'side_squared':str(cost),'violation':str(m-cost),'label_strictly_eligible_mass_above_1':True},
         'whole_SW_three_atom_counter_formula':'4sqrt(2)+sqrt(5)/2-25/4',
         'whole_SW_strict_rational_lower':str(whole_lower),'strict_excess_above_half':str(whole_lower-F(1,2)),
         'scaled_local_three_atom_gamma_formula':'4sqrt(2)/9+sqrt(5)/54-209/324',
         'scaled_local_three_atom_gamma_rational_upper':str(clipped_upper),
         'scaled_local_three_atom_gamma_below_m_over_2':True,
         'scope':'all closed unit caps for full S3 embedding pass, nonunit global cap fails; clipped formula belongs only to scaled local three atoms, not full S3; no complete halo union rerun',
         'receipt':{'python':sys.version,'platform':platform.platform(),'input_sha256':hashlib.sha256(raw).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'wall_seconds':time.perf_counter()-start}}
    args.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k]for k in ['all_passed','cell_count','maximum_unit_square_mass','whole_SW_strict_rational_lower','strict_excess_above_half','scaled_local_three_atom_gamma_rational_upper']}))


if __name__=='__main__':main()
