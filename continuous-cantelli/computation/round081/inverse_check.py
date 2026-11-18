"""Verify finite coefficient extraction and calibrate the inverse constants."""

import json
import sys
from pathlib import Path

import mpmath as mp
import sympy as sp


def exact_hermite_check():
    x,t = sp.symbols('x t')
    beta = sp.Rational(1,7)
    b = beta*(x*x-1)
    max_n = 16
    # Exact normal moment generating expression for this polynomial audit.
    exponent = -(1+beta)*t*t/2+t*t/(2*(1-beta*t*t))
    series = sp.series(sp.exp(exponent)/sp.sqrt(1-beta*t*t)-1,t,0,max_n+3).removeO()
    def mean(poly):
        return sum(c*sp.factorial2(k-1) for (k,),c in sp.Poly(poly,x).terms() if k%2==0)
    for n in range(max_n+1):
        # The identity is checked with unnormalized Hermites to retain rational arithmetic.
        lhs=mean(b*sp.hermite_prob(n,x))
        rhs=2*sp.factorial(n)*series.coeff(t,n+2)
        for j in range(2,(n+2)//2+1):
            q=n+2-2*j
            rhs -= sp.factorial(n)*mean(b**j*sp.hermite_prob(q,x))/(2**(j-1)*sp.factorial(j)*sp.factorial(q))
        assert sp.expand(lhs-rhs)==0
    return max_n+1


def main():
    rows = exact_hermite_check()
    with mp.workdps(65):
        eta = (mp.exp(mp.mpf('0.5'))-mp.mpf('1.5'))/mp.mpf('0.5')
        lower = 1-eta*eta
        assert lower>mp.mpf('0.91')
        derivative_majorant = 4*24*mp.exp(-mp.mpf(23)/24)
        assert derivative_majorant<48
        out = {'scope':'Exact finite coefficient algebra and scalar bound checks',
               'python':sys.version,'sympy':sp.__version__,'mpmath':mp.__version__,
               'exact_hermite_identity_rows':rows,
               'eta_half':mp.nstr(eta,42),'variance_fraction_lower_bound':mp.nstr(lower,42),
               'heat_second_derivative_majorant_constant':mp.nstr(derivative_majorant,42)}
    p=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).with_name('inverse-results.json')
    p.write_text(json.dumps(out,indent=2)+'\n')
    print(rows, out['eta_half'],out['variance_fraction_lower_bound'])


if __name__=='__main__':
    main()
