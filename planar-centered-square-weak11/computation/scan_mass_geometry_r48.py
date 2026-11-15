#!/usr/bin/env python3
"""R48 local R24 mass parameters at the R43 geometry, exact labelled tail."""
import json,os,sys
from fractions import Fraction as F
from itertools import product
from pathlib import Path
from periodic_halo import sha
from scan_supercell_mass import canon,digest
from scan_mass_response import run
P=Path('/Volumes/data/CodexTemp/T2/planar-centered-square-weak11');basepath=P/'round008/joint-groups/batch2/candidate_0611.json';targetpath=P/'round043/cross-group-refine/candidate_0195.json';binary=P/'round017/binary-mass/manifest.json';oldpath=Path('/Volumes/data/CodexProjects/T2/planar-centered-square-weak11/computation/inputs/r014_final_0935.json');familypath=Path('/Volumes/data/CodexProjects/T2/planar-centered-square-weak11/computation/inputs/r024_common_family.json');out=P/'round048/mass-geometry'
def main():
 if not os.path.ismount('/Volumes/data'):raise RuntimeError('data not mounted')
 base=json.loads(basepath.read_text());target=json.loads(targetpath.read_text());old=json.loads(oldpath.read_text());bm=json.loads(binary.read_text());family=json.loads(familypath.read_text());assert sha(basepath)==bm['base_sha256'] and sha(oldpath)==bm['target_sha256']
 vectors=[[F(x) for x in v['vector']] for v in bm['selected_vectors']];assert len(vectors)==8 and all(sum(v)==0 for v in vectors)
 orig=[F(p[2]) for p in old['pts']]; center=[F(-1,8)]*2+[F(1)]*2+[F(42,5)]*4
 def vals(c):return [w+sum((k*v[i] for k,v in zip(c,vectors)),F()) for i,w in enumerate(orig)]
 assert vals(center)==[F(p[2]) for p in target['pts']]
 expected={tuple(vals([a,a,F(1),F(1),b,b,b,b])) for a,b in product((F(-1),F(1)),(F(0),F(10)))}; actual={tuple(F(p[2]) for p in c['pts']) for c in family['corners']};assert expected==actual
 upper=[F(1001,1000)*F(base['pts'][i%len(base['pts'])][2]) for i in range(len(orig))];points=[];rejected=[]
 for i,j in product(range(-4,5),repeat=2):
  a=F(-1,8)+F(i,200);b=F(42,5)+F(j,20);c=[F(i,200)]*2+[F(0)]*2+[F(j,20)]*4; masses=vals([a,a,F(1),F(1),b,b,b,b]);
  if min(masses)<=0 or any(w>u for w,u in zip(masses,upper)) or sum(masses)!=sum(orig):rejected.append({'a':str(a),'b':str(b)});continue
  points.append(c)
 prior=P/'round043/cross-group-refine';pm=json.loads((prior/'manifest.json').read_text());rr={r['id']:r for r in json.loads((prior/'results.json').read_text())};cache={}
 for t in pm['tasks']:
  r=rr[t['id']]
  if 'error' not in r and t['input_sha256']==r['input_sha256'] and digest(t['input'])==t['input_sha256']:cache[t['input_sha256']]=(t,r,sha(prior/'results.json'))
 out.mkdir(parents=True,exist_ok=False);selection={'grid':'a=-1/8+i/200,b=42/5+j/20; i,j=-4..4','center':['-1/8','42/5'],'coefficients':'deltas from center: (da,da,0,0,db,db,db,db)','base_mass_model':str(oldpath),'base_mass_model_sha256':sha(oldpath),'source_family_sha256':sha(familypath),'source_binary_sha256':sha(binary),'driver_sha256':sha(__file__),'family_corners_exactly_reconstructed':True,'positive_and_kappa_cap_rejected':rejected,'kappa_cap':'1001/1000','no_fixed_geometry_envelope_inference':True}
 run(base,target,basepath,targetpath,out,10,0,points,vectors,cache,selection)
if __name__=='__main__':main()
