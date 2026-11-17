"""Fixed-angle numerical reconstruction of Gibbs Definition 9, not a coverage proof."""
from math import radians,pi,sqrt,atan2
import json,sys
from reconstruct_s_primitives import *
from reconstruct_ah import solve_2x2_rows
sigma=1.5494
n,m,o,d=labels(sigma);_,_,p=primitives(sigma)
R=solve_2x2_rows(m[2],n[1],.5,0)
S=min(circle_line_intersections(R,n[0]),key=lambda x:norm(sub(x,o['C1'])))
T=min(circle_circle_intersections(S,d['B3']),key=lambda x:norm(sub(x,p['J'])))
U=min(circle_line_intersections(d['C3'],unit(7*pi/6)),key=lambda x:norm(sub(x,p['J'])))
z=sub(U,S);d2=dot(z,z);departure=dot(sub(T,S),sub(U,T));assert d2>1
V=min([add(S,add(mul(1/d2,z),mul(k*sqrt(d2-1)/d2,(-z[1],z[0])))) for k in (-1,1)],key=lambda x:norm(sub(x,T)))
angle=atan2(cross(sub(T,S),sub(V,S)),dot(sub(T,S),sub(V,S)))
if departure>=0:green=cross(T,U);mode='straight'
else:green=arc_green_integral(T,V,S)+cross(V,U);mode='arc-plus-tangent'
green+=arc_green_integral(U,p['J'],d['C3'])+arc_green_integral(p['J'],T,d['B3'])
EH=abs(green/2)
c=__import__('math').cos(radians(sigma));P=sqrt(3)/2-2*sqrt(3)*(c/sqrt(3)-.5)**2/(4*c*c-3)
AH=.00008151836891801715773478188412529156
assert abs(norm(sub(S,R))-1)<1e-12 and abs(norm(sub(T,S))-1)<1e-12
assert abs(norm(sub(T,d['B3']))-1)<1e-12 and abs(norm(sub(U,d['C3']))-1)<1e-12
assert abs(dot(unit(7*pi/6),U)-.5)<1e-12
print(json.dumps({'python':sys.version,'sigma_degrees':sigma,'R':R,'S':S,'T':T,'U':U,'V':V,'J':p['J'],'U_distance_S':sqrt(d2),'departure_dot':departure,'path_mode':mode,'arc_angle':angle,'EH_area_numerical':EH,'full_area_numerical':P-sum(p['areas'].values())-AH-EH,'scope':'floating geometric area reconstruction only; not an independently proved universal-cover bound'},indent=2))
