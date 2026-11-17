"""Independent 80-digit algebraic endpoint diagnostic; not interval arithmetic."""
import json,sys
import mpmath as m
m.mp.dps=80
half=m.mpf(1)/2
sigma=m.mpf('1.5494')*m.pi/180
v=lambda a:m.matrix([m.cos(a),m.sin(a)])
n3=m.matrix([-1,0]);n4=v(4*m.pi/3);n5=v(5*m.pi/3)
def meet(a,b):return m.lu_solve(m.matrix([[a[0],a[1]],[b[0],b[1]]]),m.matrix([half,half]))
F3=meet(n3,v(5*m.pi/6+sigma));E3=meet(n4,v(7*m.pi/6+sigma))
foot=F3-(m.fdot(n5,F3)-half)*n5
perp=m.matrix([-n5[1],n5[0]])
r=m.sqrt(1-(m.fdot(n5,F3)-half)**2)
G=max([foot+r*perp,foot-r*perp],key=lambda a:a[0])
theta=m.atan2((G-F3)[1],(G-F3)[0])+m.pi/6
M=m.matrix([-m.mpf(1)/4,-m.sqrt(3)/4])
def L(s):
 x=m.cos(m.pi/6+s)-half
 return m.matrix([x,(half*x-half)/(m.sqrt(3)/2)])
def X(l,n):
 q=n-l;d=m.norm(q);h=m.sqrt(1-d*d/4)
 return (l+n)/2-h*m.matrix([-q[1],q[0]])/d
C=L(theta);Y=X(C,E3);Z=X(L(0),M)
fmt=lambda z:m.nstr(z,70)
result={'python':sys.version,'mpmath_version':m.__version__,'decimal_precision':80,'sigma_degrees':'1.5494','theta_radians':fmt(theta),'Y_distance_G':fmt(m.norm(Y-G)),'Y_distance_Ltheta':fmt(m.norm(Y-C)),'Z_distance_Ltheta':fmt(m.norm(Z-C)),'Lminus_G_residual':fmt(m.norm(L(-theta)-G)),'scope':'high precision numerical diagnostic, not directed intervals or proof of universal coverage'}
assert abs(m.norm(Y-C)-1)<m.mpf('1e-70')
assert m.norm(L(-theta)-G)<m.mpf('1e-70')
assert abs(m.norm(Y-G)-1)>m.mpf('.002')
print(json.dumps(result,indent=2))
