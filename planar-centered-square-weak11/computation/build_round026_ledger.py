#!/usr/bin/env python3
"""Check archived R26 identities, run counts, rational certificates and R24 comparison."""
from pathlib import Path
from fractions import Fraction as F
import gzip,hashlib,json,math
p=Path(__file__).resolve().parent; r=p/'results/round026'; read=lambda q:json.loads(q.read_text()); sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest()
def zipped(name):
 with gzip.open(r/(name+'.gz'),'rb') as f:return json.load(f)
m=zipped('manifest.json'); rows=zipped('results.json');s=read(r/'summary.json');d=read(r/'geometry_diagnostics.json');a=read(r/'primary_old_R10.json');b=read(r/'independent_old_R10.json');old=read(p/'inputs/r025_final_0110.json');new=read(p/'inputs/r026_final_0360.json');t=read(p/'results/round024/branch_tree.json')
assert len(rows)==len(m['tasks'])==s['points']==s['completed']==793 and s['reused']==1 and s['new_calls']==792 and not s['errors']
assert a['input_sha256']==b['input_sha256']==sha(p/'inputs/r026_final_0360.json')
for x,y in [('inner_area','area_lower'),('outer_area','area_upper'),('inner_ratio','weak_ratio_lower'),('outer_ratio','weak_ratio_upper')]:assert F(a[x])==F(b[y]['fraction'])
assert a['code_sha256']==sha(p/'periodic_halo.py') and b['script_sha256']==sha(p/'independent_periodic_halo.py')
assert sha(p/'scan_quadrant_groups.py')==m['scan_sha256'] and sha(p/'inputs/r025_final_0110.json')==m['target_sha256']
changes=[(i,tuple(F(x)-F(y) for x,y in zip(n,o))) for i,(n,o) in enumerate(zip(new['pts'],old['pts'])) if n!=o]
assert changes==[(79,(F(328734960629,50000000000000000),F(0),F(0))),(80,(F(328734960629,50000000000000000),F(0),F(0)))]
px,py,lam=map(F,(new['Px'],new['Py'],new['lam']));mass=sum(F(v[2]) for v in new['pts']);R=F(10);tail=mass/(px*py)*(1+px/(2*R))*(1+py/(2*R));assert tail==F(a['tail_bound'])<lam and math.ceil(R/min(px,py))+1==3
l=F(a['inner_ratio']); gap=l-F(t['upper']);oldgap=l-F(read(p/'results/round025/primary_final.json')['outer_ratio']);assert gap>0 and oldgap>0
assert sha(p/'results/round024/branch_tree.json')=='90b6ad00bff70c4346a91063f7e89588bde8e6e6e1c95e13a9b6699a69781108'
ledger={'scope':'file identity, selected input, exact arithmetic and counts; no full-search geometry replay','input_sha256':sha(p/'inputs/r026_final_0360.json'),'primary_sha256':sha(r/'primary_old_R10.json'),'independent_sha256':sha(r/'independent_old_R10.json'),'search_rows':793,'search_new_calls':792,'search_reused':1,'additional_diagnostic_baseline_calls':1,'additional_diagnostic_seconds':d['baseline']['elapsed_seconds'],'wall_seconds':s['wall_seconds'],'child_cpu_seconds':s['children_user_cpu_seconds']+s['children_system_cpu_seconds'],'lower':str(l),'upper':a['outer_ratio'],'gain_over_r25_upper':str(oldgap),'gain_over_r24_family_upper':str(gap),'tail_bound':str(tail),'q':3,'changed_label_indices_zero_based':[79,80]}
(r/'integration_checks.json').write_text(json.dumps(ledger,indent=2)+'\n');print(json.dumps(ledger,indent=2))
