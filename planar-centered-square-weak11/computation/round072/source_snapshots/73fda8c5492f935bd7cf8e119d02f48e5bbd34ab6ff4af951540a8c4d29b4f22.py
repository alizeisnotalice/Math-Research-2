#!/usr/bin/env python3
"""Exact integer coverage of two actual inputs on the R72 aligned partition."""
import hashlib,importlib.util,json,os,sys,time
from pathlib import Path
from fractions import Fraction as F
import numpy as np
t=time.monotonic();p=Path('computation/round072/aligned_cell_area_mip.py');sp=importlib.util.spec_from_file_location('ac',p);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);g=m.geometry()
src=Path('computation/results/round071-joint/one-column-result.json');d=json.loads(src.read_text());lookup={tuple(F(a)*1000 for a in p):F(w) for p,w in zip(d['primal']['positions'],d['primal']['weights'])}
weights=np.array([int(lookup.get(p,F(0))*10000) for p in m.O],dtype=np.int64)
assert all(F(int(a),10000)==lookup.get(p,F(0)) for a,p in zip(weights,m.O))
cases=[]
for name,q,approx in [('s25-only',np.zeros(145,dtype=np.int64),1.373993652089341),('r071-positive33',weights,1.4419810899674625)]:
 cm=g['C'].astype(np.int64)@q;sel=[]
 for k,(u,v,s,t) in enumerate(g['cells']):
  if not np.all(cm*100<=g['resnum'][k]):continue
  ds=[max(p[0]-u,v-p[0],p[1]-s,t-p[1]) for p in m.O+m.S]
  if any(sum(int(q[j])*100 for j,a in enumerate(ds[:145]) if a<=r)+40000*sum(a<=r for a in ds[145:])>=4*r*r+1 for r in set(ds) if r<500):sel.append(k)
 area=sum(F((g['cells'][i][1]-g['cells'][i][0])*(g['cells'][i][3]-g['cells'][i][2]),m.D) for i in sel)
 original=Path('/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-sol-kernel-sources/planar-centered-square-weak11/computation/results/round072-actual')/(name+'.json')
 cases.append({'name':name,'selected_cells':sel,'whole_cell_coverage':str(area),'whole_cell_coverage_float':float(area),'actual_deletion_approx_from_geometry_agent':approx,'coverage_loss_approx':approx-float(area),'exact_complete_area_source':str(original),'exact_complete_area_source_sha256':hashlib.sha256(original.read_bytes()).hexdigest(),'residual_closure_integer_checks':len(g['cells'])*len(g['rows'])})
out={'python':sys.version,'numpy':np.__version__,'pid':os.getpid(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'model_source_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'r071_primal_source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'scope':'Exact integer full-prefix/complete residual-delta coverage on 225 aligned cells; actual complete areas supplied and independently audited by geometry agents, pins retained; no cell score is relabelled actual deletion area.','cases':cases,'seconds':time.monotonic()-t}
Path('computation/results/round072-cells/partition-baseline-coverage.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'pid':os.getpid(),'seconds':out['seconds'],'cases':[{k:r[k] for k in ['name','whole_cell_coverage','coverage_loss_approx']} for r in cases]},indent=2))
