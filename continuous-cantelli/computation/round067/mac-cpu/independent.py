"""Independent direct scaled-input quadrature and exact quadratic Bessel checks."""
import json
from pathlib import Path
import mpmath as mp

mp.mp.dps=53
rows=[]
for a in map(mp.mpf,['0.4','0.8','1.5','2']):
    for k,j in [(2,1),(2,12),(2,256),(16,1),(64,6),(512,1)]:
        c=mp.mpf('.875')
        n=k*j
        L=(c*n*n/2)**(1/a)
        eta=j/L
        top=min(mp.pi*L,mp.power(240,1/a))
        edges=sorted(set([mp.mpf(0),top]+[mp.power(x,1/a) for x in (1,3,7,15,30,60,120,240) if mp.power(x,1/a)<top]))
        def shape(t):
            if t==0: return mp.mpf(1)
            theta=t/(2*L)
            z=t**a*(mp.sin(theta)/theta)**a
            return mp.exp(-z)*mp.cos(eta*t)
        val=mp.quad(shape,edges)/mp.pi
        limit=mp.gamma(1/a)/(mp.pi*a)
        info=dict(alpha=str(a),periods=k,harmonic=j,L=str(L),eta=str(eta),scaled_coefficient=str(val))
        if a==2:
            limit*=mp.exp(-eta**2/4)
            exact=L*mp.exp(-c*n*n)*mp.besseli(j,c*n*n)
            info['scaled_bessel']=str(exact)
            info['absolute_bessel_discrepancy']=str(abs(val-exact))
            assert abs(val-exact)<mp.mpf('1e-48')
        info['scaled_limit']=str(limit)
        info['absolute_limit_difference']=str(abs(val-limit))
        rows.append(info)
# A direct full-circle integral for a frequency forbidden by the period.
k,n,a=3,8,mp.mpf('.8')
phase=lambda x:mp.exp(mp.j*n*x-n*n*(2*abs(mp.sin(k*x/2)))**a/2)
forbidden=mp.quad(phase,[mp.pi*i/k for i in range(2*k+1)])/(2*mp.pi)
assert abs(forbidden)<mp.mpf('1e-46')
out=dict(digits=mp.mp.dps,mpmath=mp.__version__,rows=rows,forbidden=dict(k=k,n=n,direct=str(forbidden)))
Path('independent.json').write_text(json.dumps(out,indent=2)+"\n")
print("max Bessel disagreement",max(mp.mpf(r['absolute_bessel_discrepancy']) for r in rows if 'absolute_bessel_discrepancy' in r))
print("forbidden harmonic",forbidden)
