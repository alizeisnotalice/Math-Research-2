#!/usr/bin/env python3
"""Add explicit equality-empty and variable/cell metadata to final-cut MAT."""
import hashlib,json,os,sys
from pathlib import Path
import numpy as np
import scipy
from scipy.io import loadmat,savemat
from scipy import sparse
p=Path('computation/results/round072-cells/aligned225.mat');data={k:v for k,v in loadmat(p).items() if not k.startswith('__')}
model=json.loads(Path('computation/results/round072-cells/aligned225-model-exact.json').read_text());n=len(model['outside_positions_units']);nc=len(model['cells_units']);nh=len(model['labels']);nv=data['A'].shape[1]
assert nv==n+nc+nh
data.update(Aeq=sparse.csc_matrix((0,nv)),beq=np.zeros((0,1)),outside_mass_variable_count=n,cell_binary_variable_count=nc,label_binary_variable_count=nh,units_denominator=model['units_denominator'],rhs_denominator=model['rhs_denominator'],outside_positions_units=np.array(model['outside_positions_units']),S_positions_units=np.array(model['S_positions_units']),cells_units=np.array(model['cells_units']),variable_order='outsideweights[145],cell_z[225],full_label_h[340]',scipy_objective=float((data['f'].T@data['scipy_x'])[0,0]),scipy_certified_cell_area=1.383636,original_mat_sha256=hashlib.sha256(p.read_bytes()).hexdigest())
out=p.with_name('aligned225-matlab.mat');savemat(out,data,do_compression=True)
receipt={'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'pid':os.getpid(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'input_mat_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'output_mat_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'rows':data['A'].shape[0],'columns':nv,'outside_mass_variables':n,'cell_binary_variables':nc,'label_binary_variables':nh,'intcon_count':data['intcon'].size,'scipy_objective':data['scipy_objective'],'scipy_certified_cell_area':data['scipy_certified_cell_area'],'note':'Matrix/schema export only; no MILP rerun or change to exact primal. MATLAB final-cut solve is a numerical cross-check, not completion of ungenerated residual rows.'}
p.with_name('matlab-export-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
