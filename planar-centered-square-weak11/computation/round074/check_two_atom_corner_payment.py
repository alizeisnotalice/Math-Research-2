#!/usr/bin/env python3
"""Exact checks for a true two-atom corner payment; no candidate search.

The continuum inequalities are proved in the accompanying note. Finite tests
include all stated rational regression cases and one admissible n=3 input.
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


def length(a,b,c,d):return max(F(0),min(b,d)-max(a,c))


def complete_closures(points,masses):
    axes=[sorted({p[j]for p in points})for j in (0,1)];masks=set()
    for i,a in enumerate(axes[0]):
        for b in axes[0][i:]:
            for j,c in enumerate(axes[1]):
                for d in axes[1][j:]:
                    masks.add(sum(1<<k for k,p in enumerate(points)if a<=p[0]<=b and c<=p[1]<=d))
    rows=[]
    for mask in sorted(masks-{0}):
        ids=[i for i in range(len(points))if mask>>i&1]
        span=max(max(points[i][j]for i in ids)-min(points[i][j]for i in ids)for j in (0,1))
        mass=sum((masses[i]for i in ids),F(0));cap=max(F(1),span*span)
        rows.append({'mask':mask,'indices':ids,'mass':str(mass),'span_squared':str(span*span),
                     'eligible':mass>span*span,'square_cap_slack':str(cap-mass)})
    return rows


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();start=time.perf_counter();raw=a.input.read_bytes();d=json.loads(raw)
    counts1=counts2=0;min_slack=None;sharp=None
    lo,hi=d['centre_offset_integer_min'],d['centre_offset_integer_max'];den=d['centre_offset_denominator']
    for u in map(F,d['interval_lengths']):
        offsets=[F(k,den)*u for k in range(lo,hi+1)]
        for ratio in map(F,d['inner_radius_ratios']):
            v=ratio*u;delta=u-v
            for q in offsets:
                Lu=length(-u,0,q-u,q+u);Lv=length(-u,0,q-v,q+v)
                assert 0<=Lv<=Lu<=u and Lu-Lv<=delta
                counts1+=1
            for qx,qy in product(offsets,repeat=2):
                ux=length(-u,0,qx-u,qx+u);uy=length(-u,0,qy-u,qy+u)
                vx=length(-u,0,qx-v,qx+v);vy=length(-u,0,qy-v,qy+v)
                actual=ux*uy-vx*vy;bound=u*u-v*v;slack=bound-actual
                assert 0<=actual<=bound
                counts2+=1
                if min_slack is None or slack<min_slack:min_slack=slack
                if slack==0 and 0<v<u:
                    sharp={'u':str(u),'v':str(v),'q':[str(qx),str(qy)],'annular_increment':str(actual),'bound':str(bound)}
    g=d['global_grid_example'];n=g['n'];assert n==3
    e=F(1,n*n);h=F(n+1,n*n);pts=[(i*h,j*h)for i in range(n)for j in range(n)]
    pts.append(tuple(map(F,g['outside_point'])));mass=[e]*(n*n)+[F(g['outside_mass'])]
    rows=complete_closures(pts,mass);assert all(F(r['square_cap_slack'])>=0 for r in rows)
    eligible=[r for r in rows if r['eligible']];L=max(F(r['mass'])for r in eligible)
    winners=[r['mask']for r in eligible if F(r['mass'])==L]
    assert L==1 and winners==[(1<<(n*n))-1]
    parent=list(range(len(pts)))
    def find(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    for r in eligible:
        for i in r['indices'][1:]:parent[find(i)]=find(r['indices'][0])
    assert len({find(i)for i in range(len(pts))})==1
    span=max(max(p[j]for p in pts)-min(p[j]for p in pts)for j in (0,1))
    gap=span*span-F(35,33)*sum(mass,F(0));assert gap>0
    P_hi=F(1,2)-e;corner_internal=[i for i,p in enumerate(pts[:n*n])if -1<p[0]<P_hi and -1<p[1]<P_hi]
    assert corner_internal==[0] and e/2==F(1,18)
    out={'schema':'r074-two-atom-corner-payment-result-v1','all_passed':True,
         'one_dimensional_regression_cases':counts1,'two_dimensional_regression_cases':counts2,
         'minimum_annulus_bound_slack':str(min_slack),'annulus_sharp_equality_example':sharp,
         'global_n3_example':{'source_count':len(pts),'closures':rows,'maximum_mass':str(L),'maximum_masks':winners,
                              'connected':True,'nonterminal_gap':str(gap),'effective_corner_internal_indices':corner_internal,
                              'theorem_corner_deletion_upper':'1/18','complete_halo_area_computed':False,
                              'positions':[[str(v)for v in p]for p in pts],'masses':list(map(str,mass))},
         'scope':'analytic all-position two-atom true residual payment; finite rational regressions are checks, not the general proof; no bound for multiple effective outside atoms claimed',
         'receipt':{'python':sys.version,'platform':platform.platform(),'input_sha256':hashlib.sha256(raw).hexdigest(),
                    'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'wall_seconds':time.perf_counter()-start}}
    a.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k]for k in ['all_passed','one_dimensional_regression_cases','two_dimensional_regression_cases','minimum_annulus_bound_slack']}))


if __name__=='__main__':main()
