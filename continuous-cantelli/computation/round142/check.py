from fractions import Fraction as Q
import json
m,H,a,r=map(Q,['19/20','21/20','9/20','9/10'])
margins=[m-2*a,m-2*a+a*a,(r+a-1)**2-(H-m),r*r-(H-m)]
assert all(x>0 for x in margins) and 0<1-a<r<m
print(json.dumps(dict(coarse_center=str(margins[0]),center=str(margins[1]),middle=str(margins[2]),far=str(margins[3]))))
