#!/usr/bin/env python3
"""Verify R43 archived identities, selected changes, exact bounds and counts."""
from pathlib import Path
from fractions import Fraction as F
import gzip,hashlib,json,math
p=Path(__file__).resolve().parent;r=p/'results/round043';read=lambda q:json.loads(q.read_text());sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest()
def zipped(n):
 with gzip.open(r/(n+'.gz'),'rb') as f:return json.load(f)
m=zipped('manifest.json');rows=zipped('results.json');s=read(r/'summary.json');a=read(r/'primary_old_R10.json');b=read(r/'independent_old_R10.json');old=read(p/'inputs/r041_final_0036.json');new=read(p/'inputs/r043_final_0195.json');previous=read(p/'results/round041/primary_old_R10.json')
assert len(rows)==len(m['tasks'])==s['points']==s['completed']==241 and s['new_calls']==238 and s['reused']==3 and not s['errors']
assert a['input_sha256']==b['input_sha256']==sha(p/'inputs/r043_final_0195.json') and a['code_sha256']==sha(p/'periodic_halo.py') and b['script_sha256']==sha(p/'independent_periodic_halo.py')
assert sha(p/'scan_cross_group_refine_r42.py')==m['script_sha256'] and sha(p/'inputs/r041_final_0036.json')==m['target_sha256']
for x,y in [('inner_area','area_lower'),('outer_area','area_upper'),('inner_ratio','weak_ratio_lower'),('outer_ratio','weak_ratio_upper')]:assert F(a[x])==F(b[y]['fraction'])
changes=[(i,tuple(F(x)-F(y) for x,y in zip(n,o))) for i,(n,o) in enumerate(zip(new['pts'],old['pts'])) if n!=o];indices=[26,27,77,78,128,129,179,180];assert changes==[(i,(F(-3616084566919,1440000000000000000),F(0),F(0))) for i in indices]
px,py,lam=map(F,(new['Px'],new['Py'],new['lam']));mass=sum(F(v[2]) for v in new['pts']);R=F(10);tail=mass/(px*py)*(1+px/(2*R))*(1+py/(2*R));assert tail==F(a['tail_bound'])<lam and math.ceil(R/min(px,py))+1==3
lower=F(a['inner_ratio']);gap=lower-F(previous['outer_ratio']);assert gap>0
out={'scope':'file identity, selected input, exact arithmetic and counts; no full-search geometry replay','input_sha256':sha(p/'inputs/r043_final_0195.json'),'primary_sha256':sha(r/'primary_old_R10.json'),'independent_sha256':sha(r/'independent_old_R10.json'),'search_rows':241,'search_new_calls':238,'reused':3,'wall_seconds':s['wall_seconds'],'child_cpu_seconds':s['child_user_cpu_seconds']+s['child_system_cpu_seconds'],'lower':str(lower),'upper':a['outer_ratio'],'gain_over_r41_upper':str(gap),'tail_bound':str(tail),'q':3,'changed_label_indices_zero_based':indices}
(r/'integration_checks.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
