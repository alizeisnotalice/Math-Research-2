#!/usr/bin/env python3
"""R39 independent tiny sweep tests and saved active-box scalar replay.
Never enumerates candidate source boxes or imports the slab enumerator.
"""
import argparse
from decimal import Decimal, localcontext
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import importlib
import json
from math import isqrt,lcm
from pathlib import Path
import platform
import subprocess
import sys
import time

sys.dont_write_bytecode=True
parser=argparse.ArgumentParser(description=__doc__)
for name in ['operations','slab','slope','family','temp','output']:
    parser.add_argument('--'+name,required=True,type=Path)
parser.add_argument('--toy-slabs',type=Path)
args=parser.parse_args()
for key,value in vars(args).items():
    if value is not None:setattr(args,key,value.resolve())
args.temp.mkdir(parents=True,exist_ok=True)
start=time.process_time()
sys.path.insert(0,str(args.operations))
prod=importlib.import_module('r039_endpoint_sweep')
primitive=importlib.import_module('r039_interval_primitives')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def conv(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {str(k):conv(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [conv(v) for v in x]
    return x
@lru_cache(None)
def root(q,bits):
    D=1<<bits
    # Alternate exact floor identity: floor(sqrt(n*d*D^2))/d,
    # rather than production floor(sqrt(floor(n*D^2/d))).
    k=isqrt(q.numerator*q.denominator*D*D)//q.denominator
    assert k*k*q.denominator<=q.numerator*D*D<(k+1)*(k+1)*q.denominator
    lo=F(k,D);return lo,(lo if lo*lo==q else F(k+1,D))
def canon(e):
    a,v,sgn,q=e;a,v,q=F(a),F(v),F(q)
    assert q>=0 and sgn in [-1,0,1]
    if not sgn or not q:return a,v,0,F(0)
    n,d=isqrt(q.numerator),isqrt(q.denominator)
    if n*n==q.numerator and d*d==q.denominator:return a+sgn*F(n,d),v,0,F(0)
    return a,v,sgn,q
def enclosure(e,bits):
    a,v,sgn,q=e
    if not sgn:return a,a
    lo,hi=root(q,bits)
    return (a+lo,a+hi) if sgn>0 else (a-hi,a-lo)
def merge(intervals):
    z=sorted((a,b) for a,b in intervals if a<b)
    if not z:return 0
    out=0;a,b=z[0]
    for c,d in z[1:]:
        if c>b:out+=b-a;a,b=c,d
        else:b=max(b,d)
    return out+b-a

tests={};root_cases=0
for q in [F(0),F(1,9),F(2),F(4),F(10**50+1,17),F(1,10**100)]:
    for bits in [0,4,80]:
        assert primitive.sqrt_interval(q,bits)==root(q,bits);root_cases+=1
tests['independent_root_enclosures']=root_cases
assert prod.expr(0,1,1,F(4,9))==(F(2,3),F(1),0,F(0))
assert prod.compare(prod.expr(1,1,-1,2),prod.expr(1,1,-1,2))==0
rejected=[]
def expect_reject(name,fn):
    try:fn()
    except ValueError:rejected.append(name);return
    raise AssertionError(name+' was accepted')
expect_reject('nonpersistent_irrational_equality',lambda:prod.compare(prod.expr(1,0,1,2),prod.expr(1,1,1,2)))
a=prod.expr(0,0,1,2);b=prod.expr(F(1,1<<100),0,1,2)
expect_reject('unresolved_low_bits',lambda:prod.compare(a,b,80))
assert prod.compare(a,b,160)==-1
E=prod.expr
r=prod.slope_certificate([(E('-1/10',1),E('19/10',1),E(0),E(1))],0,10,0,2)
assert r['endpoint_radius']==F(1,20) and r['beta_lower']==r['beta_upper']==1
tests['raw_clip_crossing_radius']=r['endpoint_radius']
r=prod.slope_certificate([(E(1),E(3,1),E(0),E(2)),(E(1),E(3,1),E(1),E(3))],0,10,0,5)
assert r['beta_lower']==r['beta_upper']==3
tests['persistent_group_overlap_beta']=r['beta_lower']
r=prod.slope_certificate([(E(1),E(2,1),E(0),E(1)),(E(2,1),E(4),E(0),E(3))],0,10,0,5)
assert r['beta_lower']==r['beta_upper']==-2
tests['same_class_enter_leave_negative_beta']=r['beta_lower']
expect_reject('nonconstant_y_velocity',lambda:prod.slope_certificate([(E(1),E(2),E(0,1),E(1))],0,5,0,5))
expect_reject('empty_x_clip',lambda:prod.slope_certificate([(E(-3),E(-2),E(0),E(1))],0,5,0,5))
thin=[(E(1),E(2,1),E(0,0,1,2),E(F(1,1<<100),0,1,2))]
expect_reject('thin_y_at_insufficient_bits',lambda:prod.slope_certificate(thin,0,5,0,5,bits=80))
r=prod.slope_certificate(thin,0,5,0,5,bits=160)
assert 0<r['beta_lower']<=F(1,1<<100)<=r['beta_upper']
tests['thin_y_high_bits_positive_enclosure']=(r['beta_lower'],r['beta_upper'])
tests['fail_closed_cases']=rejected
# Replay only tiny saved slabs through the actual current join CLI.
toy=[]
for label,expected in [('positive',F(2)),('negative',F(-1)),('dead',F(0))]:
    if args.toy_slabs is not None:
        input_path=args.toy_slabs/(label+'.json')
    else:
        old=json.loads((args.operations/('r039-toy-'+label+'.json')).read_text())
        input_path=Path(old['argv'][1])
    dest=args.temp/('current-join-'+label+'.json')
    run=subprocess.run([sys.executable,'-B',str(args.operations/'r039_join_slope.py'),str(input_path),'--output',str(dest)],capture_output=True,text=True)
    assert run.returncode==0,run.stderr
    d=json.loads(dest.read_text());assert F(d['beta_lower'])==F(d['beta_upper'])==expected
    assert 'conditional' in d['scope'] and d['script_sha256']==sha(args.operations/'r039_join_slope.py')
    slab=json.loads(input_path.read_text());end=d['sweep']['endpoint_radius']
    assert F(d['delta'])==min(F(slab['H']),F(slab['delta_empty']),F(end) if end is not None else F(slab['H']))
    toy.append({'case':label,'beta':str(expected),'delta':d['delta'],'input_sha256':sha(input_path),'new_output_sha256':sha(dest)})
tests['current_join_tiny_replays']=toy

# Saved real active-box replay. This does not enumerate or validate rejected boxes.
slab_raw=args.slab.read_bytes();slab=json.loads(slab_raw)
slope_raw=args.slope.read_bytes();slope=json.loads(slope_raw)
assert slope['slab_sha256']==hashlib.sha256(slab_raw).hexdigest()
assert slope['slab_code_sha256']==slab['code_sha256']
assert slope['script_sha256']==sha(args.operations/'r039_join_slope.py')
for name,h in slope['dependency_sha256'].items():assert sha(args.operations/name)==h
assert slope['s0']==slab['s0'] and slope['tail']==slab['tail']
assert 'conditional' in slope['scope']
bits=slope['sweep']['bits'];D=1<<bits
xl,xr,yl,yr=map(F,slab['domain']);left,right=canon([xl,0,0,0]),canon([xr,0,0,0])
boxes=[tuple(canon(e) for e in box) for box in slab['boxes']]
assert len(boxes)==slab['counts']['active_boxes']==len(slab['box_records'])==slope['sweep']['box_count']
allx={left,right}
for l,r,y0,y1 in boxes:
    assert y0[1]==y1[1]==0;allx.update((l,r))
rows=sorted((enclosure(e,bits)[0],enclosure(e,bits)[1],e) for e in allx)
assert len(rows)==slope['sweep']['endpoint_count']
rank={row[2]:i for i,row in enumerate(rows)}
endpoint_radius=None;limiter=None
for a,b in zip(rows,rows[1:]):
    gap=b[0]-a[1];assert gap>0
    dv=abs(b[2][1]-a[2][1])
    if dv:
        v=gap/(2*dv)
        if endpoint_radius is None or v<endpoint_radius:endpoint_radius,limiter=v,(a[2],b[2])
assert endpoint_radius==F(slope['sweep']['endpoint_radius'])
assert limiter==tuple(canon(e) for e in slope['sweep']['limiting_pair'])
events=[[] for _ in rows];ybounds=[];yscale=1
for i,(l,r,y0,y1) in enumerate(boxes):
    lo=max(rank[l],rank[left]);hi=min(rank[r],rank[right]);assert lo<hi
    a,b=enclosure(y0,bits);c,d=enclosure(y1,bits)
    inner=max(yl,b),min(yr,c);outer=max(yl,a),min(yr,d)
    assert inner[0]<inner[1] and outer[0]<=inner[0]<inner[1]<=outer[1]
    ybounds.append((inner,outer));events[lo].append((1,i));events[hi].append((-1,i))
    for z in (*inner,*outer):yscale=lcm(yscale,z.denominator)
yints=[tuple(tuple(int(z*yscale) for z in interval) for interval in pair) for pair in ybounds]
active=set();beta_lo=beta_hi=F(0);gap_checks=0;max_active=0
for j,row in enumerate(rows):
    for sign,i in events[j]:
        if sign>0:assert i not in active;active.add(i)
        else:assert i in active;active.remove(i)
    if j==len(rows)-1:break
    dv=rows[j+1][2][1]-row[2][1]
    if dv:
        gap_checks+=1;max_active=max(max_active,len(active))
        lo=F(merge(yints[i][0] for i in active),yscale)
        hi=F(merge(yints[i][1] for i in active),yscale)
        assert lo<=hi
        beta_lo+=dv*(lo if dv>0 else hi);beta_hi+=dv*(hi if dv>0 else lo)
assert not active
assert beta_lo==F(slope['beta_lower'])==F(slope['sweep']['beta_lower'])
assert beta_hi==F(slope['beta_upper'])==F(slope['sweep']['beta_upper'])
assert 0<beta_lo<=beta_hi
H,delta_empty,delta=map(F,[slab['H'],slab['delta_empty'],slope['delta']])
assert delta==min(H,delta_empty,endpoint_radius)>0
family=json.loads(args.family.read_text());base_px=F(family['base']['Px'])
assert family['copy_factors']==[2,2] and 2*base_px==xr-xl
mass,lam=F(slab['tail']['mass']),F(slab['tail']['lambda'])
assert F(slab['tail']['bound'])<lam and mass>0
normalized={'definition':'tau=s/base_Px; J=lambda*A/mass','base_Px':base_px,
            's0_tau':F(slab['s0'])/base_px,'H_tau':H/base_px,'delta_tau':delta/base_px,
            'dJ_ds_bounds':(lam*beta_lo/mass,lam*beta_hi/mass),
            'dJ_dtau_bounds':(base_px*lam*beta_lo/mass,base_px*lam*beta_hi/mass)}
with localcontext() as ctx:
    ctx.prec=30
    normalized_decimals={key:[str(Decimal(v.numerator)/Decimal(v.denominator)) for v in vals]
                         for key,vals in [('dJ_dtau',normalized['dJ_dtau_bounds']),('beta',(beta_lo,beta_hi)),('delta_tau',(normalized['delta_tau'],))]}
result={'status':'passed','scope':'join/sweep source plus independent active-endpoint/y-union replay; conditional on independent slab completeness audit',
        'not_recomputed':['candidate bbox enumeration','rejected-box loop','periodic halo areas'],
        'tests':tests,'real':{'slab_sha256':hashlib.sha256(slab_raw).hexdigest(),'slope_sha256':hashlib.sha256(slope_raw).hexdigest(),
            'slab_code_sha256':slab['code_sha256'],'boxes':len(boxes),'raw_endpoint_classes':len(rows),
            'all_raw_endpoint_intervals_strictly_ordered':True,'all_y_inner_intervals_positive':True,
            'endpoint_radius':endpoint_radius,'limiting_pair':limiter,'H':H,'delta_empty':delta_empty,'delta':delta,
            'y_union_method':'direct integer interval sort/merge, no production UnionTree','tested_nonzero_velocity_gaps':gap_checks,
            'max_active_at_tested_gap':max_active,'beta_lower':beta_lo,'beta_upper':beta_hi,
            'beta_width':beta_hi-beta_lo,'normalization':normalized,'decimals_display_only':normalized_decimals},
        'source_sha256':{name:sha(args.operations/name) for name in ['r039_join_slope.py','r039_endpoint_sweep.py','r039_interval_primitives.py']},
        'normalization_family_sha256':sha(args.family),'python':sys.version,'platform':platform.platform(),
        'main_cpu_seconds':time.process_time()-start,'script_sha256':sha(Path(__file__)),
        'actual_command':[sys.executable,str(Path(__file__).resolve())]+[z for k,v in vars(args).items() if v is not None for z in ('--'+k,str(v))]}
args.output.write_text(json.dumps(conv(result),indent=2)+'\n')
print(json.dumps(conv({k:result[k] for k in ['status','scope','real','main_cpu_seconds','script_sha256']}),indent=2))
