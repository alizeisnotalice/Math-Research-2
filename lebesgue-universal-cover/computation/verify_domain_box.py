"""Independent standard-library verification of the outer-domain certificate."""
from pathlib import Path
from fractions import Fraction as F
import json,sys
from rational_trig import PI,scale,trig
p=Path(__file__).with_name('outer_domain_certificate.json');r=json.loads(p.read_text());D=r['disk_polygon_denominator'];points=[tuple(x) for x in r['disk_polygon_points']]
assert type(D)==int and all(type(v)==int for q in points for v in q)
assert D==10**15 and len(points)==512 and all(x*x+y*y<=F(D*D,4) for x,y in points)
def cross(a,b,c):return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
def hull(P):
 P=sorted(set(P));parts=[]
 for seq in (P,list(reversed(P))):
  L=[]
  for v in seq:
   while len(L)>1 and cross(L[-2],L[-1],v)<=0:L.pop()
   L.append(v)
  parts+=L[:-1]
 return parts
mins=[]
for row,n in zip(r['regular_checks'],(3,5,7,9,11,19)):
 assert row['body']=='regular_'+str(n);w=F(row['translation_coordinate_halfwidth']);assert w==F(21,100)
 s=F(row['point_radius_lower']);co=trig(scale(PI,F(1,2*n)),'cos');cn=trig(scale(PI,F(1,n)),'cos');R=1/(2*co[1]);g=R*R+w*w+2*R*w*cn[0];assert g>s*s
 H=hull(points+[(int(s*D),0)]);A=F(sum(a[0]*b[1]-a[1]*b[0] for a,b in zip(H,H[1:]+H[:1])),2*D*D)
 assert A==F(row['inscribed_polygon_area']) and A>F(419,500);mins.append(A)
row=r['fourier_check'];w=F(row['translation_coordinate_halfwidth']);assert w==F(6,25)
a=F(1,2)-F(1,80)-F(19,1200);assert a==F(283,600)==F(row['support_lower'])
s=a+w;assert s==F(row['point_radius_lower'])
H=hull(points+[(int(s*D),0)]);A=F(sum(x[0]*y[1]-x[1]*y[0] for x,y in zip(H,H[1:]+H[:1])),2*D*D)
assert A==F(row['inscribed_polygon_area']) and A>F(419,500);mins.append(A)
assert len(r['regular_checks'])==6 and r['tau']=='419/500'
print(json.dumps({'python':sys.version,'verified_regular_types':6,'verified_fourier':True,'minimum_outer_area_lower':str(min(mins)),'tau':'419/500','status':'passed','scope':'excludes outside normalized translation box only, not its interior'},indent=2))
