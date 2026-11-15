#!/usr/bin/env python3
"""Independent split-at-source 256-bit integration of the finite eight-atom test set."""
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
import datetime
import hashlib
import json
import math
import os
import platform
import subprocess
import sys
import time
import flint
from flint import arb,ctx
ctx.prec=256
AXIS=[(Q(1),Q(0)),(Q(-1),Q(0)),(Q(0),Q(1)),(Q(0),Q(-1))]
DIAG=[(Q(x,2),Q(y,2)) for x in [-1,1] for y in [-1,1]]


def fs(x):return str(x.numerator)+'/'+str(x.denominator)

def A(x):return arb(x.numerator)/x.denominator

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def endpoint(x,lower=True):
    return Q(str((x.lower() if lower else x.upper()).fmpq()))


def record(x):
    lo=endpoint(x);hi=endpoint(x,False)
    assert lo<=hi
    return {'ball':str(x),'lower_rational':fs(lo),'upper_rational':fs(hi)}


@lru_cache(None)
def primitive(t):
    assert t>=0
    if t==0:return arb(0)
    value=A(t)
    return value*value.log()-value


def cut_interval(lo,hi,source):
    points=[lo]+([source] if lo<source<hi else [])+[hi]
    return [(a,b) for a,b in zip(points,points[1:])]


def integral_u(box,source):
    # Independent implementation: physically split at source lines; on each
    # subrectangle both absolute coordinates are monotone. No signed H primitive.
    xl,xh,yl,yh=box
    sx,sy=source
    total=arb(0)
    for a,b in cut_interval(xl,xh,sx):
        ul,uh=sorted([abs(a-sx),abs(b-sx)])
        for c,d in cut_interval(yl,yh,sy):
            vl,vh=sorted([abs(c-sy),abs(d-sy)])
            total+=primitive(uh+vh)-primitive(ul+vh)-primitive(uh+vl)+primitive(ul+vl)
    return total


def reflect_interval(lo,hi,sign):
    return (lo,hi) if sign==1 else (-hi,-lo)


def expand_geometry(entries):
    expanded=[];central=0;arms=0
    for entry in entries:
        box=tuple(map(Q,entry['box']));xl,xh,yl,yh=box
        assert 0<=xl<xh<=32 and 0<=yl<yh<=1
        mult=entry['multiplicity']
        if xh<=1:
            assert mult==4
            central+=1; swapped=[False]
        else:
            assert xl>=1 and mult==8
            arms+=1;swapped=[False,True]
        orbit=[]
        for swap in swapped:
            xx=(yl,yh) if swap else (xl,xh)
            yy=(xl,xh) if swap else (yl,yh)
            for xs in [-1,1]:
                for ys in [-1,1]:
                    orbit.append((*reflect_interval(*xx,xs),*reflect_interval(*yy,ys)))
        assert len(orbit)==mult
        expanded.extend(orbit)
    scale=math.lcm(*(x.denominator for box in expanded for x in box))
    integer_boxes=[tuple(int(x*scale) for x in box) for box in expanded]
    assert all(Q(v,scale)==x for ib,box in zip(integer_boxes,expanded) for v,x in zip(ib,box))
    overlap=0
    for i,a in enumerate(integer_boxes):
        for b in integer_boxes[i+1:]:
            if max(a[0],b[0])<min(a[1],b[1]) and max(a[2],b[2])<min(a[3],b[3]):
                overlap+=1
    assert overlap==0
    original_area=sum(Q(e['multiplicity'])*(Q(e['box'][1])-Q(e['box'][0]))*(Q(e['box'][3])-Q(e['box'][2])) for e in entries)
    actual_area=sum((b-a)*(d-c) for a,b,c,d in expanded)
    assert original_area==actual_area
    return {'input_rectangles':len(entries),'central_rectangles':central,'arm_rectangles':arms,
            'expanded_rectangles':len(expanded),'expanded_pairs_checked':len(expanded)*(len(expanded)-1)//2,
            'positive_area_overlaps':overlap,'exact_integer_scaling':scale,
            'expanded_set_area':fs(actual_area),
            'reason':'Central rectangles reflected in four quadrants; arm rectangles reflected and swapped into eight sectors. All interiors are disjoint; boundary overlaps are null.'}


def main():
    assert os.path.ismount('/Volumes/data'), 'data volume not mounted'
    t0=time.monotonic()
    base=Path(__file__).resolve().parent.parent
    inp=base/'computation/inputs/cauchy-eight-atom-mass-floor/rectangles.json'
    data=json.loads(inp.read_text())
    entries=data['rectangles'];assert len(entries)==187
    geometry=expand_geometry(entries)
    va=arb(0);vd=arb(0);pieces=[]
    for i,entry in enumerate(entries):
        box=tuple(map(Q,entry['box']));mult=entry['multiplicity']
        origin=integral_u(box,(Q(0),Q(0)))
        ia=origin-sum((integral_u(box,z) for z in AXIS),arb(0))/4
        id_=origin-sum((integral_u(box,z) for z in DIAG),arb(0))/4
        va+=mult*ia;vd+=mult*id_
        pieces.append({'index':i,'box':[fs(x) for x in box],'multiplicity':mult,
                       'axis_integral_before_multiplicity':record(ia),
                       'diagonal_integral_before_multiplicity':record(id_)})
    threshold=Q(40,19)
    assert endpoint(va)>threshold and endpoint(vd)>threshold
    lower_mass=min(endpoint(va),endpoint(vd))
    # Exact geometric factor is >19/10; this finite mass witness forces cost >4.
    assert lower_mass*Q(19,10)>4
    gitbase=subprocess.run(['git','-C',str(base.parent),'rev-parse','HEAD'],
                           capture_output=True,text=True,check=True).stdout.strip()
    out={'status':'PASS_all_theta_mass_gt_40_over_19_independent_256bit_split_source',
         'timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
         'python':sys.version,'executable':sys.executable,'flint':flint.__version__,
         'precision_bits':ctx.prec,'workers':1,'platform':platform.platform(),'git_base':gitbase,
         'script_sha256':sha(Path(__file__).resolve()),'input_sha256':sha(inp),
         'author_script_sha256':data['author_script_sha256'],
         'author_certificate_sha256':data['author_certificate_sha256'],
         'original_interval_values_used':False,'original_signed_primitive_imported':False,
         'integral_algorithm':'Split every rectangle at source x/y lines, map each part to a positive distance box, use F(uh+vh)-F(ul+vh)-F(uh+vl)+F(ul+vl), F(t)=tlogt-t,F0=0.',
         'geometry':geometry,'mass_threshold':fs(threshold),
         'axis_integral_on_S':record(va),'diagonal_integral_on_S':record(vd),
         'uniform_mass_lower_rational':fs(lower_mass),
         'uniform_cost_lower_using_factor_19_over_10':fs(lower_mass*Q(19,10)),
         'pieces':pieces,'elapsed_seconds':time.monotonic()-t0,
         'scope':'Finite-set lower witness for all theta in [0,1], not a truncation of support or exact full mass. Geometry and normalization-factor proof are in companion theorem. No numerical theta interpolation, no large 2D integration, no MATLAB or LP.'}
    dest=base/'computation/results/cauchy-eight-atom-mass-floor-independent/checks-256bit.json'
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ['status','geometry','axis_integral_on_S','diagonal_integral_on_S','uniform_cost_lower_using_factor_19_over_10','elapsed_seconds']},indent=2))

if __name__=='__main__':main()
