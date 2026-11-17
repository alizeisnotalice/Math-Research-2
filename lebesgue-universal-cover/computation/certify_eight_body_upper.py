"""Rational certificate for one eight-body placement; NOT a universal cover."""
from fractions import Fraction as F
from math import factorial, isqrt
from functools import lru_cache
import json,sys,hashlib
from pathlib import Path
S=10**9; T=10**18
K=int(sys.argv[1]) if len(sys.argv)>1 else 128
TARGET=F(sys.argv[2]) if len(sys.argv)>2 else F(2097,2500)
assert K>=64
PARAMS=['-0.001596257929712958', '-0.002919442393794826', '0.05121464896977835', '0.02464167646388727', '-0.005976067545118195', '0.7189979342752278', '0.004516203660704413', '-0.012407545755070916', '0.34972813881990994', '-0.0005916505031997194', '-0.0007999287826682254', '0.1136138969315887', '-0.0011450052586001469', '0.006069465985436147', '0.23180815743002847', '-0.0006532317393185777', '-0.0021340749559243728', '2.0876915299317496', '0.007169937979170267', '-0.013467590404175424']
def add(a,b):return (a[0]+b[0],a[1]+b[1])
def neg(a):return (-a[1],-a[0])
def mul(a,b):
    c=[x*y for x in a for y in b];return min(c),max(c)
def scale(a,b):return mul(a,(b,b))
def inv(a):
    assert a[0]>0;return 1/a[1],1/a[0]
def arctan_inv(q,N):
    v=sum((F((-1)**j,(2*j+1)*q**(2*j+1)) for j in range(N)),F(0))
    r=F(1,(2*N+1)*q**(2*N+1))
    return (v,v+r) if N%2==0 else (v-r,v)
# Machin identity: tan(2a)=5/12<1, hence 0<4a<pi/2;
# tan(4a-b)=1 with b=atan(1/239) and 0<4a-b<pi/2.
t2=2*F(1,5)/(1-F(1,5)**2);assert t2==F(5,12)<1
t4=2*t2/(1-t2*t2);assert t4==F(120,119)>F(1,239)
assert (t4-F(1,239))/(1+t4*F(1,239))==1
pi_machin=add(scale(arctan_inv(5,32),F(16)),scale(arctan_inv(239,8),F(-4)))
PI=(F(3141592653589793238,T),F(3141592653589793239,T))
assert PI[0]<pi_machin[0]<=pi_machin[1]<PI[1]
assert PI[1]<F(22,7)
@lru_cache(None)
def trig(a,kind):
    mid=(a[0]+a[1])/2
    # Rational centre rounded to a fixed denominator, with exact displacement.
    x=F((mid*T).__floor__(),T); rad=max(abs(a[0]-x),abs(a[1]-x))
    assert abs(x)<12
    term=x if kind=='sin' else F(1);value=term
    for j in range(1,41):
        denominator=(2*j)*(2*j+1) if kind=='sin' else (2*j-1)*(2*j)
        term*= -x*x/denominator;value+=term
    # Polynomial includes degree 81 (sin) or 80 (cos); degree-81 Taylor
    # polynomial for cos has its zero odd term. Derivatives are bounded by 1.
    remainder=abs(x)**82/F(factorial(82))+rad
    lo=value-remainder;hi=value+remainder
    # Outward rational rounding keeps downstream arithmetic small.
    return F((lo*T).__floor__(),T),F((hi*T).__ceil__(),T)
def point(x,y):
    out=[]
    for v in (x,y):
        q=(((v[0]+v[1])/2)*S+F(1,2)).__floor__()
        assert F(q-1,S)<=v[0]<=v[1]<=F(q+1,S)
        out.append(q)
    return tuple(out)
pts=[]
for j in range(8*K):
    a=scale(PI,F(j,4*K));pts.append(point(scale(trig(a,'cos'),F(1,2)),scale(trig(a,'sin'),F(1,2))))
p=list(map(F,PARAMS));placements=[(3,F(0),p[0],p[1])]+[(n,p[2+3*i],p[3+3*i],p[4+3*i]) for i,n in enumerate((5,7,9))]
for n,r,tx,ty in placements:
    radius=scale(inv(trig(scale(PI,F(1,2*n)),'cos')),F(1,2))
    for l in range(n):
        phi=add((r,r),scale(PI,F(2*l,n)))
        vx=add((tx,tx),mul(radius,trig(phi,'cos')));vy=add((ty,ty),mul(radius,trig(phi,'sin')))
        for j in range(K+1):
            a=add(phi,scale(PI,F(1)-F(1,2*n)+F(j,n*K)))
            pts.append(point(add(vx,trig(a,'cos')),add(vy,trig(a,'sin'))))
# Additional diameter-one regular eleven-gon: its edges need no arc padding.
r11=scale(inv(trig(scale(PI,F(1,22)),'cos')),F(1,2))
for j in range(11):
    a=add((p[11],p[11]),scale(PI,F(2*j,11)))
    pts.append(point(add((p[12],p[12]),mul(r11,trig(a,'cos'))),add((p[13],p[13]),mul(r11,trig(a,'sin')))))
# Second additional polygon P19, exact vertices with no arc padding.
r19=scale(inv(trig(scale(PI,F(1,38)),'cos')),F(1,2))
for j in range(19):
    a=add((p[14],p[14]),scale(PI,F(2*j,19)))
    pts.append(point(add((p[15],p[15]),mul(r19,trig(a,'cos'))),add((p[16],p[16]),mul(r19,trig(a,'sin')))))
# Smooth width-one Fourier body: material angle a, world angle a+rho.
assert F(49,50)**2+F(11,5)**2<F(5,2)**2
for j in range(8*K):
    def material_trig(m,kind):return trig(scale(PI,F(2*j*m,8*K)%2),kind)
    h=add((F(1,2),F(1,2)),add(scale(material_trig(3,'cos'),F(1,80)),scale(material_trig(5,'cos'),-F(19,1200))))
    hp=add(scale(material_trig(3,'sin'),-F(3,80)),scale(material_trig(5,'sin'),F(19,240)))
    world=add((p[17],p[17]),scale(PI,F(2*j,8*K)))
    co=trig(world,'cos');si=trig(world,'sin')
    x=add((p[18],p[18]),add(mul(h,co),neg(mul(hp,si))))
    y=add((p[19],p[19]),add(mul(h,si),mul(hp,co)))
    pts.append(point(x,y))
def cross(a,b,c):return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
def hull(P):
    P=sorted(set(P));parts=[]
    for seq in (P,list(reversed(P))):
        H=[]
        for a in seq:
            while len(H)>1 and cross(H[-2],H[-1],a)<=0:H.pop()
            H.append(a)
        parts.append(H[:-1])
    return parts[0]+parts[1]
H=hull(pts)
assert all(cross(a,b,p)>=0 for a,b in zip(H,H[1:]+H[:1]) for p in pts)
assert all(x*x+y*y<=4*S*S for x,y in H)
area=F(sum(a[0]*b[1]-a[1]*b[0] for a,b in zip(H,H[1:]+H[:1])),2*S*S)
# Fourier boundary has sup norm of second derivative <5/2; N=8K gives error <=5pi^2/(256 K^2), also bounding the other arcs.
# Each rational endpoint is within sqrt(2)/S < 2/S of the exact endpoint.
epsilon=5*F(22,7)**2/F(256*K*K)+F(2,S)
# Exact upward edge lengths: ceil(sqrt(integer))/S.
perimeter_upper=F(0)
for a,b in zip(H,H[1:]+H[:1]):
    d=(a[0]-b[0])**2+(a[1]-b[1])**2
    z=isqrt(d);z+=int(z*z<d)
    assert z*z>=d
    perimeter_upper+=F(z,S)
assert perimeter_upper<F(88,7)
upper=area+perimeter_upper*epsilon+F(22,7)*epsilon*epsilon
assert upper<TARGET
result={'python':sys.version,'script_sha256':globals().get('STREAMED_SOURCE_SHA256') or hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'K':K,'error_model':'Fourier sup second derivative <5/2, N=8K','coordinate_denominator':S,'parameter_decimals_interpreted_exactly':PARAMS,'sample_count':len(pts),'hull_vertices':len(H),'rational_polygon_area':str(area),'epsilon':str(epsilon),'perimeter_upper':str(perimeter_upper),'strict_upper':str(upper),'strict_upper_display_only':float(upper),'target':str(TARGET),'verified':True,'scope':'upper bound on m({D,B3,B5,B7,B9,P11,P19,Fourier}), not on Lambda','hull_grid_points':H}
print(json.dumps(result,indent=2))
