#!/usr/bin/env python3
"""Certify a raw x-endpoint crossing with overlapping y supports by rational radicals."""
import hashlib,json,sys
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
p=Path(__file__).resolve().parent;ds=[json.loads((p/n/'slab.json').read_text()) for n in ['one_axis_t0','one_axis','one_axis_t1']];keys=[(86,88,29,29),(85,88,29,29)];boxes=[]
for d in ds:
 lookup={tuple(r['bbox']):(r,b) for r,b in zip(d['box_records'],d['boxes'])};boxes.append([lookup[k] for k in keys])
A0,B0=boxes[0][0][1],boxes[0][1][1];Am,Bm=boxes[1][0][1],boxes[1][1][1];A1,B1=boxes[2][0][1],boxes[2][1][1]
assert A0[0][2]==-1 and B0[1][2]==1 and F(A0[0][1])==0 and F(B0[1][1])==1
assert A0==Am==A1 and all(F(Bm[i][3])*2==F(B0[i][3])+F(B1[i][3]) for i in range(4));assert all(F(B0[i][3])==F(B0[0][3]) for i in range(4));assert F(B1[0][3])!=F(B0[0][3])
N=1<<128;q=F(A0[0][3]);k=isqrt(q.numerator*N*N//q.denominator);l=F(k,N);h=F(k+1,N);assert l*l<=q<h*h
# A_left=a−sqrt(qA), B_right=b+sqrt(qB(t)); equality sqrt(qB)=a−b−sqrt(qA).
a,b=F(A0[0][0]),F(B0[1][0]);c0,c1=a-b-h,a-b-l;assert c0>0 and c1>0
q0,q1=F(B0[0][3]),F(B1[0][3]);dr=q1-q0;v=[(c0*c0-q0)/dr,(c1*c1-q0)/dr];lower,upper=min(v),max(v);assert 0<lower<upper<F(1,2)
xlo,xhi=F(A0[0][0])-h,F(A0[0][0])-l;dom=list(map(F,ds[0]['domain']));assert dom[0]<xlo<xhi<dom[1]
yc=F(ds[0]['y_classes'][29]);assert dom[2]<yc<dom[3] and all(F(box[2][3])>0 for row in boxes for _,box in row)
out={'scope':'Raw opposite x-endpoints of two active complete boxes cross at a unique t; y supports both contain a common interior y. Does not prove an exposed union wall or a 2D common stratum.','z_normalized':'1/1000000','boxes':[list(k) for k in keys],'t_lower':str(lower),'t_upper':str(upper),'t_lower_float':float(lower),'t_upper_float':float(upper),'x_lower':str(xlo),'x_upper':str(xhi),'common_y':str(yc),'sqrt_qA_lower':str(l),'sqrt_qA_upper':str(h),'qB_t0':str(q0),'qB_t1':str(q1),'qB_derivative':str(dr),'source_slab_sha256':[hashlib.sha256((p/n/'slab.json').read_bytes()).hexdigest() for n in ['one_axis_t0','one_axis','one_axis_t1']],'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'python':sys.version};f=p/'wall_witness.json';f.write_text(json.dumps(out,indent=2)+'\n');print('wall t',float(lower),float(upper),'in (0,1/2); output sha',hashlib.sha256(f.read_bytes()).hexdigest())
