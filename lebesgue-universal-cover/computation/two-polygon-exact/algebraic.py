"""Express triple-contact hull area in Q(u), using u's exact quadratic relation."""
import json,math
from pathlib import Path
from fractions import Fraction as F
from flint import arb,ctx
ctx.prec=192;p=Path(__file__).resolve().parent;d=json.loads((p/'micro-scan.json').read_text());tr=json.loads((p/'triple.json').read_text());H=[tuple(map(F,x)) for x in d['hullP']];Q=[tuple(map(F,x)) for x in d['hullQ']];aa,bb,cc=map(F,tr['quadratic']);U=-bb/aa;V=-cc/aa # u²=Uu+V
class K:
 def __init__(self,a=0,b=0):self.a,self.b=F(a),F(b)
 def __add__(x,y):
  if not isinstance(y,K):y=K(y)
  return K(x.a+y.a,x.b+y.b)
 __radd__=__add__
 def __neg__(x):return K(-x.a,-x.b)
 def __sub__(x,y):return x+(-y if isinstance(y,K) else K(-y))
 def __rsub__(x,y):return K(y)-x
 def __mul__(x,y):
  if not isinstance(y,K):y=K(y)
  return K(x.a*y.a+x.b*y.b*V,x.a*y.b+x.b*y.a+x.b*y.b*U)
 __rmul__=__mul__
 def inv(x):
  z=K(x.a+x.b*U,-x.b);norm=x*z;assert norm.b==0 and norm.a!=0;return K(z.a/norm.a,z.b/norm.a)
 def __truediv__(x,y):return x*(y if isinstance(y,K) else K(y)).inv()
u=K(0,1);co=(1-u*u)/(1+u*u);si=2*u/(1+u*u);rot=lambda p:(co*p[0]-si*p[1],si*p[0]+co*p[1]);dot=lambda x,y:x[0]*y[0]+x[1]*y[1];cross=lambda x,y:x[0]*y[1]-x[1]*y[0];norm=lambda P:[(b[1]-a[1],a[0]-b[0]) for a,b in zip(P,P[1:]+P[:1])];hn=norm(H);qn=norm(Q);n=[tuple(K(x) for x in v) for v in hn]+[rot(v) for v in qn];h=[tuple(K(x) for x in v) for v in H];q=[rot(v) for v in Q];u0=float(arb(tr['tan_half_root_interval']).mid());ev=lambda x:float(x.a)+float(x.b)*u0
bs=[];indices=[]
for k,nn in enumerate(n):
 hi=k if k<len(H) else max(range(len(H)),key=lambda i:ev(dot(nn,h[i])))
 qi=max(range(len(Q)),key=lambda i:ev(dot(nn,q[i]))) if k<len(H) else k-len(H)
 bs.append(dot(nn,q[qi])-dot(nn,h[hi]));indices.append([hi,qi])
ub=arb(tr['tan_half_root_interval']);val=lambda x:arb(str(x.a))+arb(str(x.b))*ub
for k in json.loads((p/'local-exact.json').read_text())['positive_hinges']:
 hi,qi=indices[k];nn=n[k]
 if k<len(H):assert all(val(dot(nn,q[qi])-dot(nn,x))>0 for j,x in enumerate(q) if j!=qi)
 else:assert all(val(dot(nn,h[hi])-dot(nn,x))>0 for j,x in enumerate(h) if j!=hi)
i,j=0,7;det=cross(n[i],n[j]);t=((-bs[i]*n[j][1]+bs[j]*n[i][1])/det,(-n[i][0]*bs[j]+n[j][0]*bs[i])/det);assert (bs[26]+dot(n[26],t)).a==0 and (bs[26]+dot(n[26],t)).b==0
ones=json.loads((p/'local-exact.json').read_text())['positive_hinges'];area=sum(F(x*y1-y*x1,2) for (x,y),(x1,y1) in zip(H,H[1:]+H[:1]));value=K(area)+sum((bs[k]+dot(n[k],t) for k in ones),K(0))/2
ub=arb(tr['tan_half_root_interval']);actual=arb(str(value.a))+arb(str(value.b))*ub;assert actual<arb(tr['area_upper']).upper();out={'minimal_area':'a + b*u','a':str(value.a),'b':str(value.b),'u_quadratic_coefficients':tr['quadratic'],'u_isolating_interval':tr['tan_half_root_interval'],'area_interval':str(actual),'active_edge_support_indices':{str(k):indices[k] for k in ones},'three_contact_edges':[0,7,26]};(p/'algebraic.json').write_text(json.dumps(out,indent=2)+'\n');print('area',actual,'positive edges',ones)
