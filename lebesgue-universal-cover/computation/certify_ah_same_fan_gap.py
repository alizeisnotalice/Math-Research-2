"""Reproduce the signed improvement for the six-direction root difference over all 75 critical parameter cells."""
import argparse,json,sys
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
from certify_ah_taller_floor_chart import geometry
from certify_ah_local_physical_support import support_bound as independent
from certify_ah_same_fan_cancellation import support_bound as coupled


def run(archive,out):
    ns,_,_=geometry(archive,[(F('.427'),F('.0068'))],far=F('.22'),negative=F('.238'),eps=F('.0073'))
    A=ns['A'];n0=ns['av'](ns['d']['n'][0]);old=new=improved=A(1);rows=[]
    terms=[(A(s*(-1)**k),ns['rotate'](n0,ns['trig'](A(h+60*k)))) for s,h in ((1,F(23)),(-1,F(36))) for k in range(3)]
    for i in range(75):
        a=F(i,7500);b=a+F(1,7500);E,v,e,b0,_=independent(ns,0,a,b);_,_,_,b1,_=coupled(ns,0,a,b)
        nominal=ns['F_at'](v,e,F(36))-ns['F_at'](v,e,F(23));u0=(nominal-E*b0(terms)).lower();u1=(nominal-E*b1(terms)).lower()
        old=min(old,u0);new=min(new,u1);improved=min(improved,(u1-u0).lower());rows.append([i,str(u0),str(u1)])
    assert old<0 and new>A('.00102896509') and improved>A('.00203985770')
    obj={'line':['.427','.0068'],'physical_parameters':['.22','.238','.0073'],'quantity':'F(36 degrees)-F(23 degrees)','reference_bin':0,'lower_bounds_before_and_after_by_subcell':rows,'former_minimum':str(old),'new_minimum':str(new),'minimum_pointwise_improvement_of_lower_endpoint':str(improved),'python':sys.version,'flint':ns['flint_version'],'bits':128,'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest()}
    out.write_text(json.dumps(obj,indent=2)+'\n')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('archive',type=Path);ap.add_argument('output',type=Path);opt=ap.parse_args();run(opt.archive.resolve(),opt.output)
