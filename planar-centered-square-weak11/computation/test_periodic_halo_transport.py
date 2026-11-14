#!/usr/bin/env python3
"""Deterministic exact checks of matching, wrapping, collisions and bounds."""
import copy
import json
import sys
from fractions import Fraction as F
from pathlib import Path
import periodic_halo as original
import periodic_halo_transport as transport

out=Path(sys.argv[1]); base={'Px':'2','Py':'2','lam':'1','pts':[['0','0','1/100'],['1/10','0','1/200']]}
def copies(nx,ny):
    p=copy.deepcopy(base); p['Px']=str(nx*F(base['Px'])); p['Py']=str(ny*F(base['Py']))
    p['pts']=[[str(F(x)+a*F(base['Px'])),str(F(y)+b*F(base['Py'])),w] for a in range(nx) for b in range(ny) for x,y,w in base['pts']]
    return p
def reject(name, data, factors, delta, labels=None):
    try: transport.calculate_transport(data,base,factors,F(delta),labels=labels)
    except ValueError as e: return {'rejected':True,'message':str(e)}
    raise AssertionError(name+' accepted')

res={}; root=original.calculate(base)
for factors in ((2,1),(1,2)):
    result=transport.calculate_transport(copies(*factors),base,factors,F(0))
    assert F(result['inner_area'])==2*F(root['inner_area'])
    assert F(result['outer_area'])==2*F(root['outer_area'])
    res['zero_displacement_'+str(factors)]={'inner_area':result['inner_area'],'outer_area':result['outer_area']}
wrapped=copies(1,1); wrapped['pts'][0][0]='199/100'
result=transport.calculate_transport(wrapped,base,(1,1),F(1,100))
assert F(result['transport_matching']['maximum_torus_distance'])==F(1,100)
res['wrapped']={'max_distance':result['transport_matching']['maximum_torus_distance']}
collision=copies(1,1); collision['pts'][1][0]='0'
result=transport.calculate_transport(collision,base,(1,1),F(1,10)); old=original.calculate(collision)
assert F(result['inner_area'])==F(old['inner_area']) and F(result['outer_area'])==F(old['outer_area'])
assert result['window_atoms_distinct'] < original.calculate(base)['window_atoms_distinct']
res['collision']={'inner_area':result['inner_area'],'outer_area':result['outer_area']}
bad=copies(1,1); bad['pts'][0][2]='1/101';res['wrong_mass']=reject('wrong_mass',bad,(1,1),0)
res['duplicate_labels']=reject('duplicate_labels',copies(1,1),(1,1),0,[(0,0,0),(0,0,0)])
res['bad_factors']=reject('bad_factors',copies(1,1),(1,0),0)
bad=copies(1,1); bad['Px']='3';res['wrong_period']=reject('wrong_period',bad,(1,1),0)
bad=copies(1,1); bad['lam']='2';res['wrong_level']=reject('wrong_level',bad,(1,1),0)
res['distance_exceeded']=reject('distance_exceeded',wrapped,(1,1),0)
res['tail_fail']=reject('tail_fail',copies(1,1),(1,1),100)
reordered=copies(1,1); reordered['pts'].reverse()
result=transport.calculate_transport(reordered,base,(1,1),0,labels=[(0,0,1),(0,0,0)])
assert F(result['inner_area'])==F(root['inner_area']) and F(result['outer_area'])==F(root['outer_area'])
res['explicit_labels']={'reordered_equal':True}
probe=json.loads((out/'probe_transport_r5.json').read_text()); p8=json.loads((out/'probe_primary_r8.json').read_text()); i8=json.loads((out/'probe_independent_r8.json').read_text())
assert F(p8['inner_area'])==F(i8['area_lower']['fraction']) and F(p8['outer_area'])==F(i8['area_upper']['fraction'])
assert F(probe['inner_area'])>=F(p8['inner_area']) and F(probe['outer_area'])<=F(p8['outer_area'])
res['probe_r5_vs_r8']={'inner_area_equal':probe['inner_area']==p8['inner_area'],'outer_area_equal':probe['outer_area']==p8['outer_area'],'r5_inner_area':probe['inner_area'],'r5_outer_area':probe['outer_area'],'r8_inner_area':p8['inner_area'],'r8_outer_area':p8['outer_area']}
(out/'tests.json').write_text(json.dumps(res,indent=2)+'\n')
print(json.dumps(res,indent=2))
