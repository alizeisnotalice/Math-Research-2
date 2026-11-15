#!/usr/bin/env python3
"""Match full R52 active bbox records and seek an exact raw endpoint order reversal."""
import hashlib,json,sys
from decimal import Decimal,localcontext
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
p=Path(__file__).resolve().parent;names=['one_axis_t0','one_axis','one_axis_t1'];ds=[json.loads((p/n/'slab.json').read_text()) for n in names];assert ds[0]['x_classes']==ds[1]['x_classes']==ds[2]['x_classes'] and ds[0]['y_classes']==ds[1]['y_classes']==ds[2]['y_classes']; maps=[{tuple(r['bbox']):(r,b) for r,b in zip(d['box_records'],d['boxes'])} for d in ds];keys=set(maps[0]);assert all(set(x)==keys for x in maps) and len(keys)==26327
for k in keys:
 a,b,c=[F(x[k][0]['mass']) for x in maps];assert 2*b==a+c
# Sort by 70-digit decimal as a candidate generator; validate every reported sign by exact dyadic root bounds.
def eval_decimal(e):
 a,_,sg,q=e;a=F(a);q=F(q);return Decimal(a.numerator)/Decimal(a.denominator)+int(sg)*(Decimal(q.numerator)/Decimal(q.denominator)).sqrt()
def iv(e,bits=128):
 a,_,sg,q=e;a=F(a);q=F(q);k=isqrt(q.numerator*(1<<(2*bits))//q.denominator);lo=F(k,1<<bits);hi=lo if lo*lo==q else F(k+1,1<<bits);return (a+lo,a+hi) if sg==1 else (a-hi,a-lo) if sg==-1 else (a,a)
def cert_order(a,b):
 al,ah=iv(a);bl,bh=iv(b);return -1 if ah<bl else 1 if bh<al else 0
records=[]
for k in keys:
 for side in (0,1):records.append((k,side,maps[0][k][1][side],maps[1][k][1][side]))
with localcontext() as c:
 c.prec=70; order0=sorted(records,key=lambda x:eval_decimal(x[2]));order1=sorted(records,key=lambda x:eval_decimal(x[3]))
ranks={(x[0],x[1]):i for i,x in enumerate(order1)};witness=None
for left,right in zip(order0,order0[1:]):
 if ranks[(left[0],left[1])]<ranks[(right[0],right[1])]:continue
 if cert_order(left[2],right[2])==-1 and cert_order(left[3],right[3])==1:
  witness={'left_bbox':list(left[0]),'left_side':left[1],'right_bbox':list(right[0]),'right_side':right[1],'t0_endpoints':[left[2],right[2]],'t_half_endpoints':[left[3],right[3]],'t0_left_interval':list(map(str,iv(left[2]))),'t0_right_interval':list(map(str,iv(right[2]))),'t_half_left_interval':list(map(str,iv(left[3]))),'t_half_right_interval':list(map(str,iv(right[3]))),'interpretation':'A raw x-endpoint order reverses by continuity; this does not establish that the wall changes the union area.'};break
out={'scope':'Matches active bbox keys and checks affine mass at t=0,1/2,1; raw endpoint reversal only, not an effective wall or two-dimensional common stratum.','active_boxes':len(keys),'all_bbox_keys_match':True,'all_bbox_masses_affine_at_three_points':True,'source_sha256':{n:hashlib.sha256((p/n/'slab.json').read_bytes()).hexdigest() for n in names},'endpoint_reversal':witness,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'python':sys.version};f=p/'box_order_comparison.json';f.write_text(json.dumps(out,indent=2)+'\n');print('matched',len(keys),'reversal',bool(witness),'output sha',hashlib.sha256(f.read_bytes()).hexdigest());print('witness',json.dumps(witness))
