"""Independent high-precision y-coordinate check through the boundary layer."""
import json
from pathlib import Path
import mpmath as mp

mp.mp.dps=57
rows=[]
for c0,a0,H0 in [('.3','-1','.27'),('1.4','2','.6'),('.8','.7','.82')]:
    c,a,H=map(mp.mpf,(c0,a0,H0))
    r=mp.mpf('.25')
    A=H*H/(2*c)
    def f(x):
        if x==0:return -2*A*mp.exp(-A)
        return mp.exp(-A)*(-mp.expm1(A*(2*x-x*x)))/x
    C=mp.quad(f,[0,.25,1])-mp.quad(lambda x:mp.exp(-A*(1-x)**2)/x,[1,2,5,mp.inf])
    base=mp.quad(lambda y:mp.exp(-y*y/(2*c)),[0,H])
    for m in (64,256,1024,4096):
        k=mp.mpf(1)/m
        mx=r**k
        edges=sorted(set([mp.mpf(0),mx]+[H*mp.exp(k*z) for z in (-60,-12,-3,0,3,12,60) if H*mp.exp(k*z)<mx]))
        def integrand(y):
            if y==0: return 1
            u=y**m
            h_u=(H/y)**m
            A_y=y*y/(2*c)
            return mp.exp(a*u-u*u/2)*(mp.exp(-A_y)-mp.exp(-A_y*(1-h_u)**2))
        estimate=mp.quad(integrand,edges)
        diff=(estimate-base)/k
        rows.append(dict(c=c0,source=a0,H=H0,k=str(k),trace=str(estimate),base=str(base),
                         finite_correction=str(diff),limit_correction=str(H*C),correction_error=str(abs(diff-H*C))))
Path('precision.json').write_text(json.dumps(dict(digits=mp.mp.dps,rows=rows),indent=2)+"\n")
print([r for r in rows if r['k']=='0.000244140625'])
