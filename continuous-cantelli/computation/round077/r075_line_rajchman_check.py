"""Check the noninteger equal-ratio obstruction using independent direct products."""
import json
import platform
import mpmath as mp

mp.mp.dps=130
depth=155
rows=[]
for p in (2,3,5):
    B=(p+mp.sqrt(p*p+4))/2
    rho=1/B
    plus=lambda t:mp.cos(mp.pi*t)
    P=mp.fprod(plus(rho**h) for h in range(1,depth+1))
    tests=[]
    for k in (4,12,32,60):
        direct=mp.fprod(abs(plus(B**h)) for h in range(-depth,k+1))
        folded=P*mp.fprod(plus(rho**h) for h in range(1,k+1))
        tests.append(dict(k=k,direct=mp.nstr(direct,22),folded_difference=mp.nstr(abs(direct-folded),6),diff_from_limit=mp.nstr(abs(direct-P*P),6)))
    rows.append(dict(p=p,contraction=mp.nstr(rho,22),limit=mp.nstr(P*P,22),checks=tests))
print(json.dumps(dict(python=platform.python_version(),mpmath=mp.__version__,precision=mp.mp.dps,product_depth=depth,rows=rows),indent=2))
