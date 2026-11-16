#!/usr/bin/env python3
"""Remove accidental maximum-block ties by one exact outside contraction.

For fixed S mass1 and an admissible input, any other eligible mass1 block
contains positive outside mass; pure proper S subsets have mass<=24/25.
Multiplying all outside masses by a factor<1 preserves every upper row,
never creates a new eligible block, and makes S the unique maximum.
The chosen factor preserves the existing 89 strict whole-cell witnesses.
"""
import hashlib,importlib.util,json,os,sys,time
from pathlib import Path
from fractions import Fraction as F
t=time.monotonic();module=Path('computation/round072/aligned_cell_area_mip.py');sp=importlib.util.spec_from_file_location('aligned',module);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);g=m.geometry()
src=Path('computation/results/round072-cells/aligned225.json');old=json.loads(src.read_text());factor=F(99999999,100000000);q=[F(a)*factor for a in old['certificate']['outside_weights']]
cm=[sum(q[j] for j in rr) for rr in g['rows']];assert all(z<=F(int(b),m.D) for z,b in zip(cm,g['capnum']))
cells=[]
for c in old['certificate']['selected_cells']:
 k=c['cell_id'];assert all(z<=F(int(b),m.D) for z,b in zip(cm,g['resnum'][k]));r=F(c['closed_full_label_radius']);mu=sum(q[j] for j in c['label_outside_indices'])+F(len(c['label_S_indices']),25);margin=mu-4*r*r;assert margin>0
 cells.append(c|{'captured_mass':str(mu),'strict_full_margin':str(margin)})
inp={'threshold':'1','background':'Lebesgue R2','geometry':'centered axis-parallel common-radius squares','S_indices':list(range(25)),'positions':[[str(F(x,1000)),str(F(y,1000))] for x,y in m.S],'masses':['1/25']*25}
for (x,y),a in zip(m.O,q):
 if a>0:inp['positions'].append([str(F(x,1000)),str(F(y,1000))]);inp['masses'].append(str(a))
out={'python':sys.version,'pid':os.getpid(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'original_candidate_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'model_source_sha256':hashlib.sha256(module.read_bytes()).hexdigest(),'contraction_factor':str(factor),'outside_weights':[str(a) for a in q],'all_caps_exact_checked':len(cm),'all_selected_cell_residual_rows_exact_checked':len(cells)*len(cm),'selected_cells':cells,'certified_area':old['certificate']['certified_area'],'smallest_strict_full_margin':str(min(F(c['strict_full_margin']) for c in cells)),'unique_S_proof':'Original input has L1; mass decrease never activates an eligible set. Every other original mass1 block contains positive outside mass, since pure proper S subsets have mass<=24/25. Each such block strictly decreases. S stays eligible mass1. Actual support and all maximum blocks are independently re-enumerated by complete halo evaluator.','seconds':time.monotonic()-t,'scope':'One purposeful exact mass contraction, same candidate positions and 89 original certified closed cells. Actual complete deletion area is evaluated separately, not inferred from continuity.'}
Path('computation/results/round072-cells/contracted-unique.json').write_text(json.dumps(out,indent=2)+'\n');p=Path('computation/results/round072-cells/contracted-unique-actual-input.json');p.write_text(json.dumps(inp,indent=2)+'\n');print(json.dumps({'pid':os.getpid(),'input':str(p),'input_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'area':out['certified_area'],'min_margin':out['smallest_strict_full_margin'],'seconds':out['seconds']},indent=2))
