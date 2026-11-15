#!/usr/bin/env python3
"""Check R48 fixed candidate, four batch counts and exact rational certificates."""
from pathlib import Path
from fractions import Fraction as F
import gzip,hashlib,json,math
p=Path(__file__).resolve().parent;r=p/'results/round048';read=lambda q:json.loads(q.read_text());sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest()
def zipped(f):
 with gzip.open(f,'rb') as z:return json.load(z)
counts={}
for sub,points,new,reuse in [('batch0',81,80,1),('batch1',81,56,25),('batch2',81,56,25),('cross',15,12,3)]:
 d=r/sub;m=zipped(d/'manifest.json.gz');rows=zipped(d/'results.json.gz');s=read(d/'summary.json');assert len(m['tasks'])==len(rows)==s['points']==s['completed']==points and s['new_calls']==new and s['reused']==reuse and not s['errors'];counts[sub]={'points':points,'new_calls':new,'reused':reuse,'wall_seconds':s['wall_seconds']}
a=read(r/'primary_old_R10.json');b=read(r/'independent_old_R10.json');old=read(p/'inputs/r043_final_0195.json');new=read(p/'inputs/r048_final_0010.json');prev=read(p/'results/round043/primary_old_R10.json')
assert a['input_sha256']==b['input_sha256']==sha(p/'inputs/r048_final_0010.json') and a['code_sha256']==sha(p/'periodic_halo.py') and b['script_sha256']==sha(p/'independent_periodic_halo.py')
for x,y in [('inner_area','area_lower'),('outer_area','area_upper'),('inner_ratio','weak_ratio_lower'),('outer_ratio','weak_ratio_upper')]:assert F(a[x])==F(b[y]['fraction'])
assert len(new['pts'])==len(old['pts'])==204 and all(F(n[j])==F(o[j]) for n,o in zip(new['pts'],old['pts']) for j in (0,1))
changes=[i for i,(n,o) in enumerate(zip(new['pts'],old['pts'])) if F(n[2])!=F(o[2])];assert len(changes)==16
px,py,lam=map(F,(new['Px'],new['Py'],new['lam']));mass=sum(F(v[2]) for v in new['pts']);assert mass==sum(F(v[2]) for v in old['pts']) and all(F(w)>0 and 0<=F(x)<px and 0<=F(y)<py for x,y,w in new['pts']);R=F(10);tail=mass/(px*py)*(1+px/(2*R))*(1+py/(2*R));assert tail==F(a['tail_bound'])<lam and math.ceil(R/min(px,py))+1==3
lower=F(a['inner_ratio']);gap=lower-F(prev['outer_ratio']);assert gap>0
out={'scope':'fixed candidate and exact count summaries; no full-search geometry replay','counts':counts,'input_sha256':sha(p/'inputs/r048_final_0010.json'),'primary_sha256':sha(r/'primary_old_R10.json'),'independent_sha256':sha(r/'independent_old_R10.json'),'lower':str(lower),'upper':a['outer_ratio'],'gain_over_r43_upper':str(gap),'tail_bound':str(tail),'q':3,'changed_mass_label_indices_zero_based':changes}
(r/'integration_checks.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
