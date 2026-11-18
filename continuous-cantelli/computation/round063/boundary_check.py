"""Verify heat-potential derivatives and boundary integrals of log-flat profiles."""

import itertools
import json
import sys
from pathlib import Path

import mpmath as mp
import sympy as sp


K_VALUES = ('0.5','1','2')
L_VALUES = ('30','100','1000','10000')


def symbolic_checks():
    s,z = sp.symbols('s z',positive=True)
    p = sp.exp(-z*z/(2*s))/sp.sqrt(2*sp.pi*s)
    primitive = 2*s*p-z*sp.erfc(z/sp.sqrt(2*s))
    formulas = (-sp.erfc(z/sp.sqrt(2*s)), 2*p, -2*z*p/s)
    assert sp.simplify(sp.diff(primitive,s)-p)==0
    for j,value in enumerate(formulas,1):
        assert sp.simplify(sp.diff(primitive,z,j)-value)==0
    return [str(f) for f in formulas]


def boundary(k_label,L_label,dps,method):
    with mp.workdps(dps):
        k,L = mp.mpf(k_label),mp.mpf(L_label)
        offset = mp.log(8)
        upper = L-offset
        cut = mp.mpf(2)
        normal = 2/mp.sqrt(2*mp.pi)
        # For v(u)=u²/log(1/u)^k, extract r² from both error integrals.
        val2_integrand = lambda t: normal*t**(k/2)*mp.exp(-t**k/2)
        exponent = 1/k+mp.mpf('0.5')
        exact2 = normal*mp.power(2,exponent)/k*mp.gammainc(exponent,upper**k/2,mp.inf)
        # Use t^k/2 coordinates to integrate the omitted-inner second derivative.
        l_edge = upper**k/2
        scaled2 = normal*mp.power(2,exponent)/k*mp.exp(-l_edge)
        integral2 = scaled2*mp.quad(lambda y: (l_edge+y)**(exponent-1)*mp.exp(-y),
                                   (0,1,4,16,64,mp.inf),method=method)

        # The exterior third derivative grows towards the moving endpoint for k<=1.
        # For k=2 it is a fixed-Gaussian-shaped integral near the low endpoint.
        if k<=1:
            length = upper-cut
            shift = upper-L-upper**k/2
            cuts3 = sorted({mp.mpf(0), *[mp.mpf(x) for x in (1,4,16,64,256) if x<length],length})
            integral3 = normal*mp.exp(shift)*mp.quad(
                lambda x: (upper-x)**(3*k/2)*mp.exp(-x+(upper**k-(upper-x)**k)/2),
                cuts3,method=method)
        else:
            cuts3 = sorted({cut,*[mp.mpf(x) for x in (4,8,16,32,64) if x<upper],upper})
            integral3 = normal*mp.exp(-L)*mp.quad(
                lambda t: t**(3*k/2)*mp.exp(t-t**k/2),cuts3,method=method)
        assert val2_integrand(upper)>0
        return integral2,exact2,integral3


def main():
    formulas = symbolic_checks()
    rows = []
    with mp.workdps(70):
        max_ref = mp.mpf(0)
        max_cross = mp.mpf(0)
        for k,L in itertools.product(K_VALUES,L_VALUES):
            primary = boundary(k,L,70,'tanh-sinh')
            check = boundary(k,L,52,'gauss-legendre')
            cross = max(abs(a-b)/abs(a) for a,b in zip(primary,check))
            ref = abs(primary[0]-primary[1])/primary[1]
            max_ref = max(max_ref,ref)
            max_cross = max(max_cross,cross)
            rows.append({'k':k,'L':L,'omitted_second_over_r2':mp.nstr(primary[0],35),
                         'third_remainder_kernel_over_r2':mp.nstr(primary[2],35),
                         'gamma_reference_relative_error':mp.nstr(ref,12),
                         'cross_method_relative_difference':mp.nstr(cross,12)})
        assert max_ref<mp.mpf('1e-60'),max_ref
        assert max_cross<mp.mpf('1e-37'),max_cross
        out = {'scope':'Symbolic derivative and positive numerical boundary integrals for explicit reference profiles',
               'python':sys.version,'sympy':sp.__version__,'mpmath':mp.__version__,
               'derivatives_for_positive_z_exact':formulas,
               'methods_digits':{'tanh-sinh':70,'gauss-legendre':52},'rows':rows,
               'max_gamma_reference_relative_error':mp.nstr(max_ref,14),
               'max_cross_method_relative_difference':mp.nstr(max_cross,14)}
    path = Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).with_name('boundary-results.json')
    path.write_text(json.dumps(out,indent=2)+'\n')
    print(len(rows), out['max_cross_method_relative_difference'])


if __name__=='__main__':
    main()
