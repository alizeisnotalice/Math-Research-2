#!/usr/bin/env python3
"""Reconstruct R48 masses and cross geometry independently in exact rational arithmetic."""
import json,pathlib,sys
from fractions import Fraction as F
from periodic_halo import sha
from scan_supercell_2x2 import canon,digest
from scan_quadrant_groups import moved
P=pathlib.Path('/Volumes/data/CodexTemp/T2/planar-centered-square-weak11');d=P/'round048/mass-geometry';base=json.loads((P/'round008/joint-groups/batch2/candidate_0611.json').read_text());target=json.loads((P/'round043/cross-group-refine/candidate_0195.json').read_text());old=json.loads(pathlib.Path('/Volumes/data/CodexProjects/T2/planar-centered-square-weak11/computation/inputs/r014_final_0935.json').read_text());v=[[F(x) for x in z['vector']] for z in json.loads((P/'round017/binary-mass/manifest.json').read_text())['selected_vectors']];orig=[F(p[2]) for p in old['pts']];center=[F(-1,8)]*2+[F(1)]*2+[F(42,5)]*4;seen=set();masscases={};maxk=F(0);minw=None
for bi in range(3):
 p=d/f'batch{bi}';m=json.loads((p/'manifest.json').read_text());rs={r['id']:r for r in json.loads((p/'results.json').read_text())}
 for t in m['tasks']:
  for coef in t['coefficients']:
   co=list(map(F,coef));assert len(co)==8 and co[0]==co[1] and co[2]==co[3]==0 and len(set(co[4:]))==1
   full=[c+x for c,x in zip(center,co)];weights=[w+sum((c*z[i] for c,z in zip(full,v)),F(0)) for i,w in enumerate(orig)];assert all(x>0 for x in weights) and sum(weights)==sum(orig)
   assert all(F(point[2])==weights[i] and point[:2]==target['pts'][i][:2] for i,point in enumerate(t['input']['pts']))
   assert digest(t['input'])==t['input_sha256'];seen.add(t['input_sha256']);minw=min(weights) if minw is None else min(minw,min(weights)); maxk=max(maxk,F(rs[t['id']]['kappa']))
  if bi==0 and t['id']=='0072':masscases['first_grid_winner']=t['input']
  if bi==2 and t['id']=='0072':masscases['final_grid_winner']=t['input']
masscases['baseline']=target;assert len(seen)==193
p=d/'cross';m=json.loads((p/'manifest.json').read_text());rs={r['id']:r for r in json.loads((p/'results.json').read_text())};directions={'baseline':([],[],F(0)),'q2_28_minus':([2],[28,29],-F(1,500000)),'q2_28_plus':([2],[28,29],F(1,500000)),'all_26_minus':([0,1,2,3],[26,27],-F(1,500000)),'all_26_plus':([0,1,2,3],[26,27],F(1,500000))}
for t in m['tasks']:
 for a in t['aliases']:
  qs,ids,step=directions[a['geometry']];x=masscases[a['mass']]
  for q in qs:x=moved(base,x,q,0,ids,step)
  assert canon(x)==canon(t['input']) and str(step)==a['shift']
  maxk=max(maxk,F(rs[t['id']]['kappa']))
assert maxk<=F(1001,1000)
out={'scope':'Exact reconstruction of all mass-grid coordinates and three by five geometry patterns; no halo-kernel rerun.','grid_unique_inputs':len(seen),'cross_aliases':sum(len(t['aliases']) for t in m['tasks']),'all_masses_positive':True,'all_masses_exactly_conserved':True,'grid_coordinates_equal_r43':True,'cross_geometries_reconstructed':True,'minimum_grid_mass':str(minw),'max_kappa':str(maxk),'script_sha256':sha(__file__),'python':sys.version};path=d/'family_audit.json';path.write_text(json.dumps(out,indent=2)+'\n');print('passed',sha(path),out['grid_unique_inputs'],out['cross_aliases'])
