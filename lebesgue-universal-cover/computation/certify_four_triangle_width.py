"""Evaluate the exact four-direction minimax expression and its attaining placement."""
import json
import platform
from flint import arb, ctx

ctx.prec = 160
pi = arb.pi()
root2 = arb(2).sqrt()
R = 1/arb(3).sqrt()
alpha = pi/24
moment = R*(root2+2*(pi/12).cos())*alpha.cos()/(2+root2)-arb(1)/2
assert arb("0.060985526796930") < moment < arb("0.060985526796932")
A = R*alpha.cos()-arb(1)/2
B = R*(pi/8).cos()-arb(1)/2
p = (3*pi/8).cos()
q = (pi/8).cos()
shift = (B-A)/(p+q)
assert A > B > 0
print(json.dumps({"python":platform.python_version(),"precision":ctx.prec,
                  "sharp_excess":str(moment),"translation_along_3pi_over_8":str(shift),
                  "equal_positive_residual":str(A+shift*p),
                  "equal_negative_residual":str(-B+shift*q)},indent=2))
