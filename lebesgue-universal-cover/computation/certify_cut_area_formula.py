"""Exact enclosure of Gibbs fixed-angle cut-area formula, conditional on coverage."""
from pathlib import Path
from sigma_parameter import SIGMA_OUT, radians
from fractions import Fraction as F
import json,sys,runpy,contextlib,io,hashlib
from rational_trig import PI,trig,scale
base=Path(__file__).parent
with contextlib.redirect_stdout(io.StringIO()):d=runpy.run_path(str(base/'certify_ah_path.py'))
I=d['I'];add=d['add'];sub=d['sub'];mul=d['mul'];dot=d['dot'];norm2=d['norm2'];perp=d['perp'];r3=d['r3']
def cross(a,b):return a[0]*b[1]-a[1]*b[0]
def unit(q, coefficient=0):
 a=radians(q, coefficient);return (I(*trig(a,'cos')),I(*trig(a,'sin')))
def meet(a,b,ra=I(F(1,2)),rb=I(F(1,2))):
 det=cross(a,b);return ((ra*b[1]-a[1]*rb)/det,(a[0]*rb-ra*b[0])/det)
def project(c,n):return sub(c,mul(dot(n,c)-F(1,2),n))
def cl(c,n,sign):
 dist=dot(n,c)-F(1,2);return add(project(c,n),mul(sign*(1-dist*dist).sqrt(),perp(n)))
def cc(a,b,sign):
 q=sub(b,a);dsq=norm2(q)
 return add(mul(I(F(1,2)),add(a,b)),mul(sign*(1-dsq/4).sqrt()/dsq.sqrt(),perp(q)))
def atan(z):
 bound=max(abs(z.lo),abs(z.hi));assert bound<F(1,2)
 val=I(0);power=z
 for j in range(32):val=val+((-1)**j)*power/(2*j+1);power=power*z*z
 err=bound**65/65;return val+I(-err,err)
def arc(a,b,c):
 u=sub(a,c);v=sub(b,c);dp=dot(u,v);assert dp.lo>0
 return cross(c,sub(b,a))+atan(cross(u,v)/dp)
def absI(a):
 assert a.lo>0 or a.hi<0
 return a if a.lo>0 else -a
n=[unit(F(j,3)) for j in range(6)];m=[unit(F(1,6)+F(j,3), 1) for j in range(6)]
A1=meet(n[1],n[2]);B3=meet(n[1],m[0]);C3=meet(n[0],m[5]);C2=meet(n[5],m[5]);E3=meet(n[4],m[3]);E2=meet(n[3],m[3]);F3=meet(n[3],m[2])
K=project(F3,m[5]);G=cl(F3,n[5],1)
Ip=project(B3,m[3]);H=project(C3,n[3]);J=cc(B3,C3,-1)
P=project(G,n[2]);Q=project(E3,n[1]);Xas=cc(G,E3,-1)
CS=absI(cross(G,C2)+cross(C2,K)+arc(K,G,F3))/2
ES=absI(cross(H,E2)+cross(E2,Ip)+arc(Ip,J,B3)+arc(J,H,C3))/2
AS=absI(cross(P,A1)+cross(A1,Q)+arc(Q,Xas,E3)+arc(Xas,P,G))/2
R=meet(m[2],n[1],I(F(1,2)),I(0));S=cl(R,n[0],-1)
T=cc(S,B3,1);U=cl(C3,unit(F(7,6)),-1)
departure=dot(sub(T,S),sub(U,T));assert departure.lo>0
EH=absI(cross(T,U)+arc(U,J,C3)+arc(J,T,B3))/2
co=I(*trig(radians(F(0), 1),'cos'))
Parea=r3/2-2*r3*(co/r3-F(1,2))*(co/r3-F(1,2))/(4*co*co-3)
saved = sys.argv
output = io.StringIO()
try:
 sys.argv = [str(base/'certify_ah_area.py'), '2048']
 with contextlib.redirect_stdout(output):
  runpy.run_path(str(base/'certify_ah_area.py'))
finally:
 sys.argv = saved
raw = output.getvalue().encode()
ah = json.loads(raw)
assert ah['cells']==2048 and ah['sigma_degrees_interval']==SIGMA_OUT
assert ah['script_sha256']==hashlib.sha256((base/'certify_ah_area.py').read_bytes()).hexdigest()
AH=I(*map(F,ah['area_interval']))
value=Parea-CS-ES-AS-AH-EH
target=F('0.8440935944')
print(json.dumps({'python':sys.version,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'AH_certificate_sha256':hashlib.sha256(raw).hexdigest(),"sigma_degrees_interval": SIGMA_OUT,'area_terms':{k:v.out() for k,v in [('P',Parea),('CS',CS),('ES',ES),('AS',AS),('AH',AH),('EH',EH)]},'EH_departure_dot_lower':str(departure.lo),'EH_points':{k:[x.out() for x in p] for k,p in [('R',R),('S',S),('T',T),('U',U),('J',J)]},'full_formula_interval':value.out(),'full_formula_display_only':[float(value.lo),float(value.hi)],'strict_target':str(target),'strict_target_passed':value.hi<target,'upper_margin':str(target-value.hi),'passed':True,'scope':'certified area formula only; identifying it with a universal convex cover still requires all geometric and covering propositions'},indent=2))
