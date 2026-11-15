#!/usr/bin/env python3
"""Reconstruct R48 four-corner evidence and R52 midpoint candidate using Fraction."""
import hashlib,json,os,platform,sys
from fractions import Fraction as F
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from scan_supercell_2x2 import canon,digest
from scan_quadrant_groups import moved
P=Path('/Volumes/data/CodexTemp/T2/planar-centered-square-weak11'); O=Path(__file__).resolve().parent/'slice';O.mkdir(exist_ok=True);assert not (O/'manifest.json').exists()
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
basepath=P/'round008/joint-groups/batch2/candidate_0611.json';source=P/'round043/cross-group-refine/candidate_0195.json';binary=P/'round017/binary-mass/manifest.json';cross=P/'round048/mass-geometry/cross';base=read(basepath);original=read(source);bm=read(binary);m=read(cross/'manifest.json');results={r['id']:r for r in read(cross/'results.json')};tasks={t['id']:t for t in m['tasks']};vectors=[[F(x) for x in t['vector']] for t in bm['selected_vectors']];velocity=[F(3,50)*(vectors[0][i]+vectors[1][i])-F(3,5)*sum((vectors[j][i] for j in range(4,8)),F()) for i in range(204)];assert sum(velocity)==0
ids={(0,0):'0000',(0,1):'0002',(1,0):'0010',(1,1):'0012'};px=F(base['Px']);mass=sum(F(p[2]) for p in original['pts']);corners=[]
for (t,z),id in ids.items():
 data={**original,'pts':[[p[0],p[1],str(F(p[2])+t*velocity[i])] for i,p in enumerate(original['pts'])]};data=moved(base,data,2,0,[28,29],F(z,500000));q=tasks[id];r=results[id];assert canon(data)==canon(q['input']) and digest(data)==q['input_sha256']==r['input_sha256'];assert all(F(p[2])>0 for p in data['pts']) and sum(F(p[2]) for p in data['pts'])==mass; assert F(r['tail_bound'])<1
 path=O/f'corner_t{t}_z{z}.json';path.write_text(json.dumps(data,indent=2)+'\n');corners.append({'t':str(t),'normalized_z':str(F(z,500000)),'task_id':id,'canonical_sha256':digest(data),'file_sha256':sha(path),'inner_ratio':r['inner_ratio'],'outer_ratio':r['outer_ratio'],'inner_area':str(F(r['inner_ratio'])*mass),'outer_area':str(F(r['outer_ratio'])*mass),'tail_bound':r['tail_bound'],'tail_radius':r['tail_radius'],'kappa':r['kappa'],'delta':r['delta']})
r00,r01,r10,r11=[results[ids[x]] for x in [(0,0),(0,1),(1,0),(1,1)]];lo=F(r11['inner_ratio'])-F(r10['outer_ratio'])-F(r01['outer_ratio'])+F(r00['inner_ratio']);hi=F(r11['outer_ratio'])-F(r10['inner_ratio'])-F(r01['inner_ratio'])+F(r00['outer_ratio']);assert lo>0 and lo<=hi
mid={**original,'pts':[[p[0],p[1],str(F(p[2])+velocity[i]/2)] for i,p in enumerate(original['pts'])]};assert min(F(p[2]) for p in mid['pts'])>0 and sum(F(p[2]) for p in mid['pts'])==mass
midpath=O/'midmass_t_half.json';midpath.write_text(json.dumps(mid,indent=2)+'\n');velpath=O/'velocity_x.json';velpath.write_text(json.dumps([int(i in (130,131)) for i in range(204)])+'\n');mwpath=O/'velocity_mass.json';mwpath.write_text(json.dumps([str(x) for x in velocity],indent=2)+'\n')
manifest={'scope':'Four finite corner enclosures; common two-parameter stratum not yet certified. Midpoint candidate is an input, not a completed slope calculation.','source_hashes':{'base':sha(basepath),'original':sha(source),'mass_vectors':sha(binary),'cross_manifest':sha(cross/'manifest.json'),'cross_results':sha(cross/'results.json'),'script':sha(__file__)},'px':str(px),'lambda':original['lam'],'mass':str(mass),'mass_derivative':'3/50*(V0+V1)-3/5*(V4+V5+V6+V7)','sum_mass_derivative':str(sum(velocity)),'geometry_velocity_nonzero_original_labels':[130,131],'z_convention':'physical z; normalized z=z/Px','corner_certificates':corners,'mixed_J_lower':str(lo),'mixed_J_upper':str(hi),'midpoint_candidate':{'mass_t':'1/2','normalized_z':'1/1000000','physical_z':str(px/F(1000000)),'midmass_sha256':sha(midpath),'coordinate_velocity_sha256':sha(velpath),'mass_velocity_sha256':sha(mwpath)},'python':sys.version,'platform':platform.platform()};(O/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print('corners',len(corners),'positive mixed difference',str(lo),str(hi));print('output',O,'manifest sha',sha(O/'manifest.json'))
if __name__=='__main__':
 if not os.path.ismount('/Volumes/data'):raise RuntimeError('data not mounted')
