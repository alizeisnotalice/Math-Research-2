#!/usr/bin/env python3
"""Verify R25 fixed-mass position grid and original-kernel certificates."""
import hashlib,json,math,shutil
from fractions import Fraction as F
from pathlib import Path
from scan_supercell_joint import make
from scan_supercell_mass import target_delta
proj=Path(__file__).resolve().parent;root=Path('/Volumes/data/CodexTemp/T2/planar-centered-square-weak11/round025/fixedmass-position');out=proj/'results/round025';base=proj/'inputs/r008_final_0611.json';previous=proj/'inputs/r023_final_1051.json';final=proj/'inputs/r025_final_0110.json'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def digest(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
source,old,new=map(read,(base,previous,final));bx,by,lam=(F(source[k]) for k in ('Px','Py','lam'));px,py=(F(new[k]) for k in ('Px','Py'))
assert (px,py,F(new['lam']))==(2*bx,2*by,lam) and len(new['pts'])==len(old['pts'])==204
mass=sum((F(p[2]) for p in new['pts']),F(0));assert mass==sum((F(p[2]) for p in old['pts']),F(0))==4*sum((F(p[2]) for p in source['pts']),F(0))
assert all(F(w)>0 and 0<=F(x)<px and 0<=F(y)<py for x,y,w in new['pts']) and all(F(a[2])==F(b[2]) for a,b in zip(old['pts'],new['pts']))
geom,_=make(source,(F(9,8),F(27,64),F(0)));assert all(F(a[j])==F(b[j]) for a,b in zip(geom['pts'],new['pts']) for j in (0,1))
changes={name:[i for i,(a,b) in enumerate(zip(old['pts'],new['pts'])) if F(a[j])!=F(b[j])] for j,name in enumerate(('x','y','weight'))}
R=F(10);tail=mass/(px*py)*(1+px/(2*R))*(1+py/(2*R));q=math.ceil(R/min(px,py))+1
p,ind=read(out/'primary_final.json'),read(out/'independent_final.json');assert tail==F(p['tail_bound'])<lam and q==3
lo,hi=F(p['inner_ratio']),F(p['outer_ratio']);oldhi=F(read(proj/'results/round023/primary_final.json')['outer_ratio']);assert lo>oldhi
assert p['input_sha256']==ind['input_sha256']==sha(final) and p['code_sha256']==sha(proj/'periodic_halo.py') and ind['script_sha256']==sha(proj/'independent_periodic_halo.py')
assert lo==F(ind['weak_ratio_lower']['fraction']) and hi==F(ind['weak_ratio_upper']['fraction'])
m,s,r=(read(root/n) for n in ('manifest.json','summary.json','results.json'));tasks=m['tasks'];completed=[json.loads(x) for x in (root/'completed.jsonl').read_text().splitlines()]
assert len(tasks)==len(r)==len(completed)==s['points']==s['completed']==289 and not s['errors']
assert sorted(completed,key=lambda x:x['id'])==r and sha(root/'results.json')==s['results_sha256'] and digest(tasks)==m['task_list_sha256']==s['manifest_task_list_sha256']
assert {x['id']:x['input_sha256'] for x in tasks}=={x['id']:x['input_sha256'] for x in r}
assert m['base_sha256']==sha(base) and m['mass_input_sha256']==sha(previous) and m['dominated_kernel_sha256']==sha(proj/'periodic_halo_dominated.py') and m['scan_sha256']==sha(proj/'scan_fixedmass_positions.py') and m['geometry_builder_sha256']==sha(proj/'scan_supercell_joint.py')
grid=set()
for item in tasks:
 assert digest(item['input'])==item['input_sha256']
 for alias in item['aliases']:
  i,j=alias['i'],alias['j'];a,b=F(19,16)+F(i,32),F(27,64)+F(j,32);assert (F(alias['a']),F(alias['b']))==(a,b);grid.add((i,j));g,_=make(source,(a,b,F(0)))
  for row,w in zip(g['pts'],old['pts']):row[2]=w[2]
  assert g==item['input'] and F(item['delta'])==target_delta(source,g)
assert grid=={(i,j) for i in range(-8,9) for j in range(-8,9)}
winner=s['top'][0];assert winner['id']==max(r,key=lambda v:F(v['inner_ratio']))['id']=='0110' and read(root/'candidate_0110.json')==new
shutil.copyfile(root/'summary.json',out/'summary.json')
ledger={'points':len(r),'new_calls':len(r),'reused':0,'wall_seconds':s['wall_seconds'],'children_user_cpu_seconds':s['children_user_cpu_seconds'],'children_system_cpu_seconds':s['children_system_cpu_seconds'],'manifest_sha256':sha(root/'manifest.json'),'task_list_sha256':m['task_list_sha256'],'results_sha256':s['results_sha256'],'scan_sha256':m['scan_sha256'],'winner':winner}
checks={'input_sha256':sha(final),'baseline_r23_sha256':sha(previous),'mass':str(mass),'Px':str(px),'Py':str(py),'level':str(lam),'positive_mass_and_cell_valid':True,'mass_four_times_r8':True,'changes_zero_based':changes,'parameters':['9/8','27/64'],'tail_radius':str(R),'tail_bound':str(tail),'q':q,'lower':str(lo),'upper':str(hi),'r23_upper':str(oldhi),'gain_over_r23_upper':str(lo-oldhi),'result_rows':len(r),'worker_evaluations':len(r),'reused_records':0,'wall_seconds':s['wall_seconds'],'child_cpu_seconds':s['children_user_cpu_seconds']+s['children_system_cpu_seconds']}
(out/'search_ledger.json').write_text(json.dumps(ledger,indent=2)+'\n');(out/'integration_checks.json').write_text(json.dumps(checks,indent=2)+'\n');print(json.dumps(checks,indent=2))
