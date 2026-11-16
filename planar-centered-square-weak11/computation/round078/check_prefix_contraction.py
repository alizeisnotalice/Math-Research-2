#!/usr/bin/env python3
"""Scalar regressions of analytic contraction constants and the specified thick-prefix cell."""
from pathlib import Path
from fractions import Fraction as F
import json,sys,os,time,platform,hashlib
start=time.perf_counter();ip=Path(sys.argv[1]);cfg=json.loads(ip.read_text());counts={'two_dimensional_comparisons':0,'one_row_squared_comparisons':0};rows=[]
for n in cfg['n_values']:
 e=F(1,n*n);h=F(n+1,n*n);K=n//2;assert h*h>e
 for k in range(2,K+1):
  assert k*e<4*(k-1)**2*h*h;counts['one_row_squared_comparisons']+=1
  for l in range(2,K+1):
   assert k*l<=4*(k-1)*(l-1)
   assert k*l*e<4*(k-1)*h*((l-1)*h+e)
   counts['two_dimensional_comparisons']+=1
 rows.append({'n':n,'K':K})
q=cfg['thick_fixture'];n=q['n'];e=F(1,n*n);h=F(n+1,n*n);a=F(q['axis_distance']);b=F(q['each_outside_mass']);B=2*b
# Scalar candidates and monotonic differences prove the all-index minima in the note.
assert n==20 and a==F(1,2) and a+9*h<1<a+10*h
V=min(1-F(10,n),(a+10*h)**2-F(11,n))
U=min(1-100*e,(a+10*h)**2-121*e)
assert V==F(1,2) and U==F(1197,1600) and U-B==F(13,1600)>0
assert h-2*h*h>0 and h*h-e>0 and h-F(1,n)>0
lo,hi=map(F,q['cell_negative_bounds']);r=F(q['common_full_radius']);t0=F(q['cell_mindistance'])
assert hi-lo==F(1,500) and a+lo==t0 and -lo<r<F(1,2)
assert max(a+hi,7*h-lo)<r and 8*h-hi>r
assert B+64*e-4*r*r==F(1079,10000)>0
assert 4*t0*t0>B+10*e and 4*t0*t0>b+F(1,4) and 4*t0*t0>B
assert -hi>F(1,2*n) and -F(1,2)<lo<hi<-e
assert e+b>a*a
nonterminal=(1-e+a)**2-F(35,33)*(1+B);assert nonterminal>0
out={'passed':True,'counts':counts,'n_rows':rows,'fixture':{'V':str(V),'U':str(U),'joint_capacity_slack':str(U-B),'whole_cell_area':str((hi-lo)**2),'common_full_radius':str(r),'full_captured_old_count':64,'full_mass':str(B+64*e),'full_strict_margin':str(B+64*e-4*r*r),'minimum_outside_capture_radius':str(t0),'threshold_lower':str(4*t0*t0),'thin_joint_mass_upper':str(B+10*e),'single_outside_mass_upper':str(b+F(1,4)),'residual_total_mass':str(B),'nonterminal_gap':str(nonterminal)},'scope':'Analytic contraction sufficient-inequality scalar regressions; exact specified receiver-cell and monotone-cap-minimum arithmetic. Not a global halo-area calculation or complete source closure enumeration.','runtime':{'python':platform.python_version(),'pid':os.getpid(),'seconds':time.perf_counter()-start,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'input_sha256':hashlib.sha256(ip.read_bytes()).hexdigest()}}
Path(sys.argv[2]).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
