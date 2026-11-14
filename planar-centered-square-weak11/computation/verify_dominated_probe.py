#!/usr/bin/env python3
"""Check saved dominated probe against original independent R10 arithmetic."""
import hashlib,json,sys
from fractions import Fraction as F
from pathlib import Path

folder=Path(sys.argv[1]); read=lambda name:json.loads((folder/name).read_text())
p=read('probe204.json'); a=read('dominated_r5.json'); b=read('primary_r10.json'); c=read('independent_r10.json'); m=a['transport_matching'];base=m['base_input']
assert a['input_sha256']==b['input_sha256']==c['input_sha256']==hashlib.sha256((folder/'probe204.json').read_bytes()).hexdigest()
ratios=[];dis=[];px,py=F(p['Px']),F(p['Py']);bx,by=F(base['Px']),F(base['Py'])
for row,label,lift in zip(p['pts'],m['labels'],m['lifts']):
    x,y,w=map(F,row); u,v,i=label; sx,sy,sw=map(F,base['pts'][i]);ix,iy=lift['period_shift'];dx,dy=map(F,lift['displacement'])
    assert (x%px+ix*px)-(sx%bx+u*bx)==dx and (y%py+iy*py)-(sy%by+v*by)==dy
    assert F(lift['base_mass'])==sw and F(lift['target_mass'])==w and F(lift['mass_ratio'])==w/sw
    ratios.append(w/sw);dis.append(max(abs(dx),abs(dy)))
assert max(ratios)==F(m['kappa']) and max(dis)==F(m['maximum_torus_distance'])<=F(m['delta'])
assert sum(F(r[2]) for r in p['pts'])==F(m['target_mass'])==F(a['mass'])
assert F(a['tail_bound'])+F(m['tail_margin'])==F(p['lam']) and F(m['tail_margin'])>0
assert F(b['inner_area'])==F(c['area_lower']['fraction']) and F(b['outer_area'])==F(c['area_upper']['fraction'])
assert F(a['inner_area'])>=F(b['inner_area']) and F(a['outer_area'])<=F(b['outer_area'])
result={'input_sha256':a['input_sha256'],'kernel_sha256':a['code_sha256'],'matching_rows':len(ratios),'kappa':m['kappa'],'maximum_distance':m['maximum_torus_distance'],'target_mass':m['target_mass'],'tail_bound':a['tail_bound'],'tail_margin':m['tail_margin'],'r5_inner_area':a['inner_area'],'r5_outer_area':a['outer_area'],'r10_inner_area':b['inner_area'],'r10_outer_area':b['outer_area'],'r5_r10_inner_equal':F(a['inner_area'])==F(b['inner_area']),'r5_r10_outer_equal':F(a['outer_area'])==F(b['outer_area']),'r5_inner_ratio':a['inner_ratio'],'r5_outer_ratio':a['outer_ratio'],'primary_seconds':b['elapsed_seconds'],'independent_seconds':c['total_seconds'],'dominated_seconds':a['elapsed_seconds']}
(folder/'probe_checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
