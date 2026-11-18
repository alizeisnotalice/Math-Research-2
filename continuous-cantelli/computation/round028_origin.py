"""Scaled positive-side integral on L/2 <= t*x <= 6L, see error bounds in proof."""
import json
import platform
import mpmath as mp
mp.mp.dps = 50
rows = []
for b in (mp.mpf(0), mp.mpf('.5'), mp.mpf(1)):
    for L in (100, 1000, 10000):
        def integrand(z):
            y = 2*L + 2*mp.sqrt(L)*z
            ell = L-mp.log(y)
            D = ell+b*mp.log(ell)
            return (mp.exp(-(y*mp.exp(-L))**2/2) * ell**b/y
                    * mp.exp(-(y-2*D)**2/(4*D)) * 2*mp.sqrt(L/(2*mp.pi)))
        lo, hi = -mp.mpf('.75')*mp.sqrt(L), 2*mp.sqrt(L)
        cuts = [lo]+[v for v in (-8,-3,0,3,8) if lo<v<hi]+[hi]
        mass = mp.quad(integrand,cuts)
        predicted = mp.power(L,b-mp.mpf('.5'))/mp.sqrt(2)
        rows.append([str(b),L,mp.nstr(mass,21),mp.nstr(mass/predicted,16)])
print(json.dumps(dict(python=platform.python_version(),mpmath=mp.__version__,
                      digits=mp.mp.dps,columns=['b','L','central_mass','ratio'],rows=rows),indent=2))
