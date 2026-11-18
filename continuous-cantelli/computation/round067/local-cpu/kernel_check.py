"""Independent symbolic derivatives and higher-precision time integration."""
import json
from pathlib import Path
import mpmath as mp
import sympy as sp

z,s=sp.symbols("z s",positive=True)
density=sp.exp(-z*z/(2*s))/sp.sqrt(2*sp.pi*s)
heat=2*s*density-z*sp.erfc(z/sp.sqrt(2*s))
expected=(-sp.erfc(z/sp.sqrt(2*s)),2*density,-2*z*density/s)
for j,formula in enumerate(expected,1):
    assert sp.simplify(sp.diff(heat,z,j)-formula)==0

mp.mp.dps=72
samples=[]
for a in map(mp.mpf,(".001",".1",".5","1","3","6","12")):
    theory=(1+a*a)*mp.erfc(a/mp.sqrt(2))/2-a*mp.exp(-a*a/2)/mp.sqrt(2*mp.pi)
    direct=mp.quad(lambda q:mp.erfc(a/mp.sqrt(2*q))/2,[0,mp.mpf(".01"),mp.mpf(".1"),1])
    error=abs(theory-direct)
    assert error<mp.mpf("1e-65")
    samples.append({"alpha":str(a),"mass":mp.nstr(direct,32),"absolute_error":mp.nstr(error,12)})
Path("kernel-check.json").write_text(json.dumps({
    "sympy":sp.__version__,"mpmath":mp.__version__,"dps":mp.mp.dps,"records":samples
},indent=2)+"\n")
