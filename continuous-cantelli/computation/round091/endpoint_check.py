"""Symbolic and numerical checks for the Stefan endpoint and fork rates."""

import itertools
import json
import sys
from pathlib import Path

import mpmath as mp
import sympy as sp


def symbolics():
    u,b,T,z,c,kappa = sp.symbols('u b T z c kappa',positive=True)
    phi = lambda x: sp.exp(-x*x/(2*T))/sp.sqrt(2*sp.pi*T)
    q = phi(b-u)*(1-sp.exp(-2*b*u/T))
    coeff = sp.simplify(sp.diff(q,u).subs(u,0))
    assert sp.simplify(coeff-2*b*phi(b)/T)==0
    A = c*z*z/2-z
    root = sp.sqrt(kappa/A)
    action = A*u+kappa/u
    assert sp.simplify(sp.diff(action,u).subs(u,root))==0
    assert sp.simplify(action.subs(u,root)**2-4*A*kappa)==0
    return str(coeff),str(root)


def reflected_cap(z_label,dps,method):
    with mp.workdps(dps):
        z=mp.mpf(z_label)
        b,T=mp.mpf('0.7'),mp.mpf('1.2')
        q0=1/mp.sqrt(2*mp.pi*T)

        def stable(r):
            u=r/z
            density=q0*mp.exp(-(b-u)**2/(2*T))*(-mp.expm1(-2*b*u/T))
            return z*mp.exp(-r)*density

        upper=2*b*z
        cuts=tuple(sorted({mp.mpf(0),*[mp.mpf(x) for x in (1,4,16,64,128,256) if x<upper],upper}))
        value=mp.quad(stable,cuts,method=method)
        limit=2*b*q0*mp.exp(-b*b/(2*T))/T
        return value,limit


def fork_rate(z_label,factor,dps,method):
    with mp.workdps(dps):
        z=mp.mpf(z_label)
        b,c=mp.mpf('0.7'),mp.mpf('1.3')
        kappa=mp.mpf(factor)*b*b/(2*c)
        A=c*z*z/2-z
        saddle=mp.sqrt(kappa/A)
        mu=mp.sqrt(A*kappa)
        mass_scaled_bessel=2*saddle*mp.besselk(1,2*mu)*mp.exp(2*mu)
        width=1/mp.sqrt(mu)
        # Dimensionless ratio w=u/saddle. Use the centered coordinate
        # s=(w-1)/width and subtract the exact action at w=1.
        upper=mp.mpf(18)
        lower=-min(upper,mp.mpf('0.99')/width)
        cuts=tuple(sorted({lower,*[mp.mpf(x) for x in (-10,-5,-2,0,2,5,10) if lower<x<upper],upper}))

        def stable(s):
            d=width*s
            return mp.exp(-mu*d*d/(1+d))

        mass_scaled_quadrature=saddle*width*mp.quad(stable,cuts,method=method)
        ref=abs(mass_scaled_bessel-mass_scaled_quadrature)/mass_scaled_bessel
        rate=b-2*mu/z+(2*mp.log(z)+mp.log(mass_scaled_bessel))/z
        return rate,b-mp.sqrt(2*c*kappa),ref


def main():
    coefficient,saddle = symbolics()
    with mp.workdps(66):
        cap_rows=[]
        cap_max=mp.mpf(0)
        for z in ('1e2','1e3','1e4','1e5'):
            actual,limit=reflected_cap(z,66,'gauss-legendre')
            check,_=reflected_cap(z,48,'tanh-sinh')
            cross=abs(check-actual)
            cap_max=max(cap_max,cross)
            cap_rows.append({'z':z,'scaled_cap':mp.nstr(actual,35),'predicted_limit':mp.nstr(limit,35),
                             'cross_method_difference':mp.nstr(cross,14)})

        fork_rows=[]
        fork_max=mp.mpf(0)
        for z,factor in itertools.product(('1e3','1e4','1e5'),('0.5','1','2')):
            value,limit,error=fork_rate(z,factor,66,'gauss-legendre')
            fork_max=max(fork_max,error)
            fork_rows.append({'z':z,'kappa_over_threshold':factor,'effective_log_rate':mp.nstr(value,30),
                              'limit':mp.nstr(limit,30),
                              'Bessel_identity_cross_error':mp.nstr(error,14)})
        assert cap_max<mp.mpf('1e-35'),cap_max
        assert fork_max<mp.mpf('1e-36'),fork_max
        result={'scope':'Exact local symbolic checks plus floating heat-cap and scalar saddle calibrations',
                'python':sys.version,'sympy':sp.__version__,'mpmath':mp.__version__,
                'exact_reflected_linear_coefficient':coefficient,
                'exact_saddle':saddle,'cap_rows':cap_rows,'fork_rows':fork_rows,
                'max_cap_cross_method_difference':mp.nstr(cap_max,15),
                'max_fork_relative_Bessel_error':mp.nstr(fork_max,15)}
    p=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).with_name('endpoint-results.json')
    p.write_text(json.dumps(result,indent=2)+'\n')
    print(len(cap_rows),len(fork_rows),result['max_cap_cross_method_difference'],result['max_fork_relative_Bessel_error'])


if __name__=='__main__':
    main()
