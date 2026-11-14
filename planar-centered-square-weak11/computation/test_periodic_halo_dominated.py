#!/usr/bin/env python3
"""Exact tests for dominated mass transport and numerical homogeneity."""
import copy, hashlib, json, platform, sys
from fractions import Fraction as F
from pathlib import Path
import periodic_halo_dominated as d
import periodic_halo_transport as t

out=Path(sys.argv[1]); base={'Px':'2','Py':'2','lam':'1','pts':[['0','0','1/100'],['1/10','0','1/100']]}
def calc(data,b=base,**kw):return d.calculate_dominated(data,b,(1,1),F(0),**kw)
def reject(name,data,b=base,**kw):
    try:calc(data,b,**kw)
    except ValueError as e:return str(e)
    raise AssertionError(name+' accepted')
def areas(x):return (F(x['inner_area']),F(x['outer_area']))

result={}; original=t.calculate_transport(base,base,(1,1),F(0)); equal=calc(base)
assert all(original[k]==equal[k] for k in original if k!='transport_matching') and equal['transport_matching']['kappa']=='1'
result['kappa_one_matches_transport']=True
move=copy.deepcopy(base);move['pts'][0][2]='1001/100000';move['pts'][1][2]='999/100000'; r=calc(move,require_mass_conservation=True)
assert r['transport_matching']['kappa']=='1001/1000' and r['transport_matching']['mass_conserved'] and F(r['mass'])==F(1,50)
result['conserved_transfer']={'kappa':r['transport_matching']['kappa'],'mass':r['mass']}
scaled_base=copy.deepcopy(base);scaled_base['lam']='2'; doubled=copy.deepcopy(scaled_base)
for row in doubled['pts']:row[2]=str(2*F(row[2]))
scaled=calc(doubled,scaled_base)
assert areas(scaled)==areas(equal) and F(scaled['inner_ratio'])==F(equal['inner_ratio']) and not scaled['transport_matching']['mass_conserved']
result['nonconserved_homogeneity']={'same_area_and_ratio':True,'kappa':scaled['transport_matching']['kappa']}
result['require_conservation_reject']=reject('nonconservation',doubled,scaled_base,require_mass_conservation=True)
zero=copy.deepcopy(base);zero['pts'][0][2]='0';r=calc(zero); assert r['transport_matching']['kappa']=='1' and F(r['mass'])==F(1,100)
result['one_zero_target_accepted']=True
allzero=copy.deepcopy(zero);allzero['pts'][1][2]='0';result['all_zero_reject']=reject('allzero',allzero)
negative=copy.deepcopy(base);negative['pts'][0][2]='-1';result['negative_reject']=reject('negative',negative)
zero_base=copy.deepcopy(base);zero_base['pts'][0][2]='0';result['zero_base_reject']=reject('zero_base',base,zero_base)
bad=copy.deepcopy(base);bad['pts'][0][2]='1000';result['tail_fail_reject']=reject('kappa-tail',bad)
bad=copy.deepcopy(base);bad['Px']='3';result['period_reject']=reject('period',bad)
result['label_reject']=reject('labels',base,labels=[(0,0,0),(0,0,0)])
wrapped=copy.deepcopy(base);wrapped['pts'][0][0]='199/100';r=d.calculate_dominated(wrapped,base,(1,1),F(1,100))
assert r['transport_matching']['lifts'][0]['period_shift']==[-1,0] and r['transport_matching']['lifts'][0]['displacement']==['-1/100','0']
result['wrap_lift_verified']=True
result.update(python=sys.version,platform=platform.platform(),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),kernel_sha256=hashlib.sha256(Path(d.__file__).read_bytes()).hexdigest())
(out/'tests.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
