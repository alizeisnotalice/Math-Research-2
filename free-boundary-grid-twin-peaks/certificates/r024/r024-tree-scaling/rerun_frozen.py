#!/usr/bin/env python3
"""Reproduce the finite tree tables with the immutable r024 source file."""
import importlib.util,json,hashlib,time
from fractions import Fraction
from pathlib import Path

HERE=Path(__file__).resolve().parent
SRC=HERE/'tree_peak_dp_r024_frozen.py'
spec=importlib.util.spec_from_file_location('frozen',SRC)
t=importlib.util.module_from_spec(spec);assert spec.loader is not None;spec.loader.exec_module(t)

def row(kind,h):
    st=time.time()
    if kind=='single': a,b,g=t.binary_counts_fast(h);n=(1<<(h+1))-1
    else:a,b,g=t.double_binary_counts_fast(h);n=2*((1<<(h+1))-1)
    r=Fraction(a,b)
    return {'height':h,'vertices':n,'A2':str(a),'Btop2':str(b),'ratio':str(r),
            'ratio_float':float(r),'pair_orbits':g,'wall_seconds':time.time()-st}

out={'source':'tree_peak_dp_r024_frozen.py','source_sha256':hashlib.sha256(SRC.read_bytes()).hexdigest(),
     'single_h1_h9':[row('single',h) for h in range(1,10)],
     'double_h1_h7':[row('double',h) for h in range(1,8)],
     'scope':'Exact finite rational recomputation. It does not assert an all-height asymptotic.'}
(HERE/'frozen-rerun.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
