#!/usr/bin/env python3
"""R76 minimal exact algebra checks; no position search or halo integration."""
from fractions import Fraction as F
import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys
import time


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();start=time.perf_counter();raw=a.input.read_bytes();d=json.loads(raw)
    m=F(d['m']);e=F(d['corner_cut']);assert m==e==F(1,9)
    pair=(1-2*e)*m/4;singleton=(F(1,6)-e)**2;total=singleton+2*pair
    assert pair==F(7,324)and singleton==F(1,324)and total==F(5,108)
    assert m/2-total==F(1,108)>0
    algebra=[]
    for t in map(F,d['u_plus_v_values']):
        delta=m/(4*t);u=(t+delta)/2;v=(t-delta)/2
        assert 0<=v<=u<=F(1,2)and u*u-v*v==m/4
        B=max(F(0),u-e)**2-max(F(0),v-e)**2
        assert 0<=B<=pair
        algebra.append({'u':str(u),'v':str(v),'clipped_pair_bound':str(B),'uniform_bound':str(pair)})
    geometry=[]
    for case in d['geometric_cases']:
        pts=[(F(0),F(0))]+[tuple(map(F,p))for p in case['outside']]
        x0=min(p[0]for p in pts);x1=max(p[0]for p in pts);y0=min(p[1]for p in pts);y1=max(p[1]for p in pts)
        delta=max((x1-x0)/2,(y1-y0)/2,x1+e,y1+e,-F(1,2)-x0,-F(1,2)-y0,F(0))
        L=m+sum(map(F,case['masses']),F(0));active=L>4*delta*delta
        assert active==case['expected_triple_active']
        geometry.append({'name':case['name'],'bbox':list(map(str,[x0,x1,y0,y1])),'critical_corner_radius':str(delta),
                         'label_mass':str(L),'strict_triple_active':active,'global_caps_checked':False})
    out={'schema':'r076-triple-null-branch-result-v1','all_passed':True,'pair_bound':str(pair),'singleton_bound':str(singleton),
         'triple_null_corner_deletion_upper':str(total),'gap_below_m_over_2':str(m/2-total),
         'six_algebra_fixtures':algebra,'three_geometry_fixtures':geometry,
         'scope':'exact constants and strict/tie geometric predicates only; no full-grid cap qualification, position search or area integration; analytic theorem in note',
         'receipt':{'python':sys.version,'platform':platform.platform(),'input_sha256':hashlib.sha256(raw).hexdigest(),
                    'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'wall_seconds':time.perf_counter()-start}}
    a.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k]for k in ['all_passed','pair_bound','singleton_bound','triple_null_corner_deletion_upper','gap_below_m_over_2']}))


if __name__=='__main__':main()
