"""Independent high-precision direct positive-core integration."""
import json
from pathlib import Path
import mpmath as mp

mp.mp.dps=60
d0=mp.log(2)/2
rows=[]
for j,q,r in [(0,0,0),(mp.mpf(1)/4,mp.mpf('0.67'),2),(-mp.mpf('.2'),mp.mpf('1.1'),-1)]:
    def D(e):
        s=mp.log(e)
        return e+s/2+d0+(j*s+q)/e+r*s*s/e**2
    for L0 in (128,512,2048,8192,32768):
        L=mp.mpf(L0)
        ell=mp.findroot(lambda e:e+mp.log(2*D(e))-L,L-mp.log(2*L))
        d=D(ell); g=mp.diff(D,ell)
        raw=mp.exp(d-ell-mp.exp(-2*ell)/2)*mp.sqrt(2*d)/(2*(d+g))
        h=2*mp.sqrt(d)
        def primal(z):
            y=2*d+h*z
            p=D(L-mp.log(y))
            return mp.exp(p-L-(y-2*p)**2/(4*p)-(y*mp.exp(-L))**2/2)
        val=h*mp.quad(primal,[-8,-4,0,4,8])/mp.sqrt(2*mp.pi)
        first=L*L*(raw-val+1/(2*L))+(j/2+mp.mpf(1)/8)*mp.log(L)+q/2-mp.mpf(1)/2+d0/4
        rows.append(dict(j=str(j),q=str(q),r=str(r),L=L0,
                         coefficient_remainder=str(first),
                         second_saddle_ratio=str(d*(2*d*(val/raw-1)-1))))
Path("precision.json").write_text(json.dumps({"digits":mp.mp.dps,"rows":rows},indent=2)+"\n")
print(rows[-1])
