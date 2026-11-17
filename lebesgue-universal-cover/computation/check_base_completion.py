"""CAS identity and exhaustive finite cap-rotation checks; analytic proof is separate."""
import sys,json,itertools
from fractions import Fraction as F
import sympy as s
x,y,a,b,p,q=s.symbols('x y a b p q',real=True)
r2=(x-a)**2+(y-b)**2
identity=s.expand(((x-p)**2+(y-q)**2)-((a-p)**2+(b-q)**2)-r2-2*((x-a)*(a-p)+(y-b)*(b-q)))
assert identity==0
u,v=s.symbols('u v',real=True)
antipode=s.expand((x-(a-u))**2+(y-(b-v))**2-r2-2*((x-a)*u+(y-b)*v)-(u*u+v*v))
assert antipode==0
cases=[]
for choices in itertools.product((0,1,2),repeat=3):
 hit={i+(3 if c==2 else 0) for i,c in enumerate(choices) if c}
 good=[k for k in range(6) if not ({(j+k)%6 for j in hit}&{3,5})]
 assert good
 cases.append({'hit':sorted(hit),'valid_rotations_in_units_pi_over_3':good})
assert len(cases)==27
sqrtlo=F(1732050807568877293,10**18);sqrthi=F(1732050807568877294,10**18)
assert sqrtlo**2<3<sqrthi**2
pal=[2-2/sqrtlo,2-2/sqrthi];assert pal[1]<F(8453,10000)
print(json.dumps({'python':sys.version,'sympy':s.__version__,'active_ball_square_expansion':str(identity),'antipode_square_expansion':str(antipode),'Pal_area_interval':[str(v) for v in pal],'cap_hit_sets_checked':27,'minimum_available_rotations':min(len(c['valid_rotations_in_units_pi_over_3']) for c in cases),'cases':cases,'passed':True,'scope':'symbolic algebra and finite cap combinatorics; compactness/separation and IVT are proved analytically'},indent=2))
