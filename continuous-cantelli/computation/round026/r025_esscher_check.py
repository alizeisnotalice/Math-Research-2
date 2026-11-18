"""Verify the shifted exponent algebra and quantify Fourier errors in the tail proof."""
import json
import platform
import mpmath as mp
import sympy as sp

x,t,u,g,eta,sigma2=sp.symbols('x t u g eta sigma2', real=True)
s=t+sp.I*u
original= -x**2/2 + s*x-s**2*(1+g)/2
shifted=-(x-t)**2/2-t**2*g/2+sp.I*u*(x-t*(1+g))+u**2*(1+g)/2
assert sp.expand(original-shifted)==0
assert sp.expand(-t**2*g/2-(eta-g)*u**2/2+((t**2-u**2)*g+eta*u**2)/2)==0
mp.mp.dps=65
rows=[]
for tk in (16,64,256,1024):
    tv=mp.mpf(tk)
    vv=80*mp.log(tv)/tv**2
    cutoff=tv/2
    # Both spectral tails occur for |u|>t/2; factor four counts two signs/laws.
    high_exact=4*mp.sqrt(mp.pi/(2*vv))*mp.erfc(mp.sqrt(vv/2)*cutoff)
    high_upper=4*mp.exp(-vv*cutoff**2/2)/(vv*cutoff)
    low_upper=tv**(-29)
    density_error=(high_upper+low_upper)/(2*mp.pi)
    assert high_exact<high_upper and density_error<mp.mpf('1e-10')
    local_floor=mp.exp(-1/(2*(1+2*vv)))/mp.sqrt(2*mp.pi*(1+2*vv))-density_error
    rows.append({'t':tk,'spectral_tail_exact':mp.nstr(high_exact,16),
                 'spectral_tail_upper':mp.nstr(high_upper,16),
                 'density_error_upper':mp.nstr(density_error,16),
                 'density_floor_on_unit_interval':mp.nstr(local_floor,16)})
print(json.dumps({'python':platform.python_version(),'sympy':sp.__version__,'mpmath':mp.__version__,
    'mpmath_digits':mp.mp.dps,'exact_shifted_exponent_checks':2,
    'rows':rows},indent=2))
