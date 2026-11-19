"""Check the small even cases and the coefficient in the analytic tail bound."""
import json
from fractions import Fraction as F
from itertools import combinations
from flint import arb, ctx

ctx.prec = 256
normals = [(F(1),F(0)),(F(21,29),F(20,29)),(F(-12,37),F(35,37)),
           (F(-1),F(0)),(F(-12,37),F(-35,37)),(F(21,29),F(-20,29))]
def det(x,y):
    return x[0]*y[1]-x[1]*y[0]
def vertex(x,y):
    z=det(x,y)
    assert z>0
    return ((y[1]-x[1])/z,(x[0]-y[0])/z)
points=[vertex(x,y) for x,y in zip(normals,normals[1:]+normals[:1])]
for x,y in normals:
    assert x*x+y*y==1
    assert max(x*a+y*b for a,b in points)==1
diam2=max((x-a)**2+(y-b)**2 for (x,y),(a,b) in combinations(points,2))
assert diam2==F(900,169) and F(169,225)-F(3,4)==F(1,900)
pi=arb.pi()
gap8=(pi/7).cos()/(pi/14).cos()-(pi/8).cos()
assert gap8 > arb(2594)/10**7
tail=F(1,2)*(1-F(10,2400))**2-F(3,8)*F(100,81)/(1-F(10,648))
assert tail==F(941599,36748800)>0
k=F(1,8)
assert (3*k*k-k**3)/3==F(23,1536)
print(json.dumps({"hexagon_squared_diameter":str(diam2),"hexagon_min_normals":[[str(x),str(y)] for x,y in normals],"octagon_improvement":str(gap8),"tail_gap":str(tail),"even_cubic_lower_pi_squared":str(F(23,1536))},indent=2))
