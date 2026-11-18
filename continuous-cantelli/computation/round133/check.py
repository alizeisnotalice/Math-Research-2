"""Compute local tilted mass using u=zx and v=(u-2B)/sqrt(2B)."""
import json
import mpmath as m
m.mp.dps=70
rows=[]
for ell in map(m.mpf,('50','200','1000','5000')):
    for beta,c in [(m.mpf('0.25'),0),(m.mpf('0.5'),m.log(2)/2)]:
        B=ell+beta*m.log(ell)+c
        L=ell+m.log(2*B)
        def kernel(t):
            h=t/m.sqrt(2*B)
            if h<=-1:return m.mpf(0)
            l=ell-m.log1p(h)
            delta=l-ell+beta*m.log(l/ell)
            return m.exp(-t*t/2+B*(1+h)**2*delta/(B+delta))/m.sqrt(2*m.pi)*m.exp(-((2*B*(1+h))*m.exp(-L))**2/2)
        v=m.quad(kernel,[max(-m.sqrt(2*B),-14),-5,0,5,14])
        mass=v*m.sqrt(2*B)*m.exp(B-L)
        rows.append([str(ell),str(beta),m.nstr(mass,22),m.nstr(mass/(m.exp(c)/m.sqrt(2)*L**(beta-m.mpf('.5'))),18),m.nstr(B*(v-1),16)])
print(json.dumps(dict(mpmath=m.__version__,digits=m.mp.dps,window=14,columns=['ell','beta','truncated_mass','leading_ratio','scaled_saddle_correction'],rows=rows),indent=2))
