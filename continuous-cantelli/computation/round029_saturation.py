"""Integrate the standardized central window for the signed correction j."""
import json
import mpmath as mp
mp.mp.dps=75
rows=[]
for j in [-1,0,1]:
    for L in map(mp.mpf,['1e3','1e8','1e16','1e24']):
        def D(ell):
            return ell+mp.log(ell)/2+mp.log(2)/2+j*mp.log(ell)/ell
        y0=mp.findroot(lambda y:y-2*D(L-mp.log(y)),2*L)
        d0=y0/2
        def weight(z):
            y=y0+2*mp.sqrt(d0)*z
            ell=L-mp.log(y)
            d=D(ell)
            return (2*mp.sqrt(d0*ell/mp.pi)/y*mp.exp(j*mp.log(ell)/ell
                     -(y-2*d)**2/(4*d)-(y*mp.exp(-L))**2/2))
        mass=mp.quad(weight,[-15,-8,-3,0,3,8,15])
        rows.append([j,str(L),mp.nstr((mass-1)*L/mp.log(L),19),j-.25])
print(json.dumps(dict(dps=mp.mp.dps,mpmath=mp.__version__,z_window=[-15,15],
                      columns=['j','L','scaled_excess','predicted_limit'],rows=rows),indent=2))
