"""Check the absorbing half-line heat kernel and its mass by reflection."""
import json, mpmath as mp
mp.mp.dps=60
rows=[]
for t in map(mp.mpf,('.25','.5','1','2')):
    def p(y):return mp.exp(-y*y/(2*t))/mp.sqrt(2*mp.pi*t)
    mass=mp.quad(lambda y:p(y)-p(y+2),[-1,0,1,mp.inf])
    expected=mp.erf(1/mp.sqrt(2*t))
    assert abs(mass-expected)<mp.mpf('1e-55')
    rows.append({'time':str(t),'survival':mp.nstr(mass,24),'difference':str(abs(mass-expected))})
print(json.dumps(dict(mpmath=mp.__version__,checks=rows),indent=2))
