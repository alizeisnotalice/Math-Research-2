"""Independent exact area/perimeter arithmetic from certificate JSON on stdin.
Does not revalidate trig enclosure or boundary coverage; use generator and proof for those.
"""
import json,sys
from fractions import Fraction as F
from math import isqrt
r=json.load(sys.stdin);v=r['hull_grid_points'];s=r['coordinate_denominator']
a=F(0);p=F(0)
for i,u in enumerate(v):
    w=v[(i+1)%len(v)];a+=F(u[0]*w[1]-u[1]*w[0],2*s*s)
    d=sum((x-y)**2 for x,y in zip(u,w));z=isqrt(d)
    if z*z<d:z+=1
    p+=F(z,s)
assert a==F(r['rational_polygon_area']) and p==F(r['perimeter_upper'])
e=F(22,7)**2/F(72*r['K']**2)+F(2,s)
u=a+p*e+F(22,7)*e*e
assert e==F(r['epsilon']) and u==F(r['strict_upper']) and u<F(r['target'])
print(json.dumps({'python':sys.version,'verified':True,'area':str(a),'perimeter_upper':str(p),'strict_upper':str(u),'scope':'independent certificate arithmetic, not independent trig generation'},indent=2))
