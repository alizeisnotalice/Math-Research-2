"""Whole-rectangle tangent tests giving injectivity of the AH coordinate map."""
from pathlib import Path
from fractions import Fraction as F
import contextlib,io,runpy,json,sys,hashlib
with contextlib.redirect_stdout(io.StringIO()):d=runpy.run_path(str(Path(__file__).with_name('certify_ah_path.py')))
I=d['I'];V=d['V'];add=d['add'];sub=d['sub'];dot=d['dot'];norm2=d['norm2'];r3=d['r3']
C,G,E,Nm=(d[k] for k in ('C','G','E','Nm'))
fL=V(F(-1,2),r3/2);fN=V(F(1,2),r3/2)
checks=[]
for ln,L in [('C',C),('G',G)]:
 for nn,N in [('Nm',Nm),('E',E)]:
  a=norm2(sub(add(L,fL),N))-1
  b=norm2(sub(add(N,fN),L))-1
  # The G,Nm equality is proved symbolically, not by an interval straddling zero.
  if (ln,nn)!=('G','Nm'):assert a.hi<0
  assert b.hi<0
  checks.append({'L':ln,'N':nn,'L_tangent_distance_minus_one':a.out(),'N_tangent_distance_minus_one':b.out(),'L_equality_analytic':(ln,nn)==('G','Nm')})
Lbox=tuple(I(min(c.lo,g.lo),max(c.hi,g.hi)) for c,g in zip(C,G))
Nbox=tuple(I(min(c.lo,g.lo),max(c.hi,g.hi)) for c,g in zip(Nm,E))
q=sub(Nbox,Lbox)
assert q[0].hi<0
assert dot(q,d['dL']).hi<0 and dot(q,fL).lo>0
qn=sub(Lbox,Nbox)
assert dot(qn,d['dN']).hi<0 and dot(qn,fN).lo>0
# Positive quadrant tests for the path. Endpoints involving Z use analytic equalities.
L0=d['L0'];M=d['M'];T=d['T'];Z=d['Z'];Y=d['Y']
assert norm2(sub(T,L0)).hi<1
assert dot(sub(Z,M),sub(T,Z)).lo>0
arc_lower_M=[];arc_upper_L0=[]
for j in range(64):
 x=T[0]*I(F(j,64),F(j+1,64));y=C[1]+(1-(x-C[0])*(x-C[0])).sqrt();X=V(x,y)
 a=norm2(sub(X,L0));b=norm2(sub(X,M))
 assert a.hi<1 and b.lo>1
 arc_upper_L0.append(a.hi);arc_lower_M.append(b.lo)
print(json.dumps({'python':sys.version,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'corners':checks,'arc_L0_squared_distance_upper':str(max(arc_upper_L0)),'arc_M_squared_distance_lower':str(min(arc_lower_M)),'segment_M_outward_derivative_half_lower':str(dot(sub(Z,M),sub(T,Z)).lo),'passed':True,'scope':'fixed-angle coordinate injectivity and positive-quadrant path placement; no area or universal-cover theorem'},indent=2))
