#!/usr/bin/env python3
"""Check R35 archived counts, fixed candidate, two kernels and rational gain."""
from pathlib import Path
from fractions import Fraction as F
import gzip,hashlib,json,math
p=Path(__file__).resolve().parent;r=p/'results/round035';read=lambda q:json.loads(q.read_text());sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest()
def zipped(f):
 with gzip.open(f,'rb') as z:return json.load(z)
counts={}
for sub,points,new,reuse in [('axis-wall-probe',401,400,1),('axis-wall-refine',251,240,11)]:
 d=r/sub;m=zipped(d/'manifest.json.gz');rows=zipped(d/'results.json.gz');s=read(d/'summary.json');assert len(m['tasks'])==len(rows)==s['points']==s['completed']==points and s['new_calls']==new and s['reused']==reuse and not s['errors'];counts[sub]={'points':points,'new_calls':new,'reused':reuse,'wall_seconds':s['wall_seconds'],'child_cpu_seconds':s['child_user_cpu_seconds']+s['child_system_cpu_seconds']}
a=read(r/'primary_old_R10.json');b=read(r/'independent_old_R10.json');old=read(p/'inputs/r029_final_0298.json');new=read(p/'inputs/r035_final_0134.json');prev=read(p/'results/round029/primary_old_R10.json')
assert a['input_sha256']==b['input_sha256']==sha(p/'inputs/r035_final_0134.json') and a['code_sha256']==sha(p/'periodic_halo.py') and b['script_sha256']==sha(p/'independent_periodic_halo.py')
for x,y in [('inner_area','area_lower'),('outer_area','area_upper'),('inner_ratio','weak_ratio_lower'),('outer_ratio','weak_ratio_upper')]:assert F(a[x])==F(b[y]['fraction'])
changes=[(i,tuple(F(x)-F(y) for x,y in zip(n,o))) for i,(n,o) in enumerate(zip(new['pts'],old['pts'])) if n!=o];assert changes==[(i,(F(2958614645661,250000000000000000),F(0),F(0))) for i in (130,131)]
px,py,lam=map(F,(new['Px'],new['Py'],new['lam']));mass=sum(F(v[2]) for v in new['pts']);R=F(10);tail=mass/(px*py)*(1+px/(2*R))*(1+py/(2*R));assert tail==F(a['tail_bound'])<lam and math.ceil(R/min(px,py))+1==3
lower=F(a['inner_ratio']);gap=lower-F(prev['outer_ratio']);assert gap>0
out={'scope':'fixed input identity, exact arithmetic and count summaries; not a full-search or secant replay','counts':counts,'input_sha256':sha(p/'inputs/r035_final_0134.json'),'primary_sha256':sha(r/'primary_old_R10.json'),'independent_sha256':sha(r/'independent_old_R10.json'),'lower':str(lower),'upper':a['outer_ratio'],'gain_over_r29_upper':str(gap),'tail_bound':str(tail),'q':3,'changed_label_indices_zero_based':[130,131]}
(r/'integration_checks.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
