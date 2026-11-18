"""Check convex heat remainder against the exact archived three-layer mixture."""
import json
import math
import platform
import mpmath as mp
from scipy.optimize import minimize_scalar
import scipy

mp.mp.dps=60
g=lambda x: mp.exp(-x*x/2)/mp.sqrt(2*mp.pi)
Phi=lambda x: mp.erfc(-x/mp.sqrt(2))/2
def call(k,V):
    z=(2-k)/mp.sqrt(V)
    return (2-k)*Phi(z)+mp.sqrt(V)*g(z)
eps=mp.mpf(3)/5*mp.exp(-mp.mpf(16)/3)
rows=[]
for k in map(mp.mpf,(-3,0,2,6)):
    sig=eps*(call(k,mp.mpf(9)/8)-call(k,mp.mpf(5)/8))/2
    def marked(s):
        variances=[mp.mpf(5)/8+s/2,mp.mpf(9)/8-s/2]
        return eps/4*sum(g((k-2)/mp.sqrt(V))/mp.sqrt(V) for V in variances)
    a=mp.quad(lambda s:(1-s)*marked(s)/2,[0,1])
    b=mp.quad(lambda s:s*marked(s)/2,[0,1])
    assert min(a,b)>0 and max(abs(sig-a),abs(sig-b))<mp.mpf('1e-58')
    rows.append(dict(k=float(k),target=mp.nstr(sig,25),maximum_error=mp.nstr(max(abs(sig-a),abs(sig-b)),8)))
kappa=2/math.sqrt(2*math.pi)*math.exp(-.5)
def q(t):
    return kappa*(1+(t-1)*math.log1p(-1/t)) if t>1 else kappa
opt=minimize_scalar(lambda t:math.log(q(t))/t,bounds=(1.00000001,4),method='bounded')
print(json.dumps(dict(python=platform.python_version(),mpmath=mp.__version__,scipy=scipy.__version__,
                     precision=60,randomness='none',model='archived R005 exact randomized three-layer law',
                     rows=rows,kappa=kappa,centered_sign_constant=kappa*math.log(2),
                     small_amplitude_optimized_t=float(opt.x),optimized_scaled_rayleigh_factor=float(-2*opt.fun)),
                 indent=2))
