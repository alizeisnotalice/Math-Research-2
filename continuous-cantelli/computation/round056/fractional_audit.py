import json
import platform

import mpmath as mp
import sympy as sp

z = sp.symbols('z', positive=True)
m=3
chi=(z-1)**8 * (2-z)**8
psi=sp.diff(chi,z,m)
for j in range(m):
    assert sp.integrate(z**j*psi,(z,1,2)) == 0
for alpha in (sp.Rational(5,2),sp.Rational(7,2)):
    target=(-1)**(m+2)*sp.prod(alpha-j for j in range(m+2))*sp.integrate(z**(alpha-m-2)*chi,(z,1,2))
    integrand=sp.integrate(z**alpha*sp.diff(psi,z,2),(z,1,2))
    assert sp.simplify(integrand-target)==0
    assert integrand != 0


def calibrate(digits):
    with mp.workdps(digits):
        phi=lambda z: mp.exp(-z*z/2)/mp.sqrt(2*mp.pi)
        p=lambda z,v: phi(z/mp.sqrt(v))/mp.sqrt(v)
        g3=lambda z,v: -(z**3/v**3-3*z/v**2)*p(z,v)
        g5=lambda z: -(z**5-10*z**3+15*z)*phi(z)
        alpha=mp.mpf('2.5')
        model_v=lambda z: (2 if z>0 else 1)*abs(z)**alpha
        cutoff=[-4,-2,-1,0,1,2,4]
        expected=phi(0)*mp.quad(lambda z:model_v(z)*g5(z)/2,cutoff)
        weak=[]
        for n in (12,32,64,96):
            h=mp.power(2,-n)
            r=h**(alpha-2)
            result=mp.quad(lambda z:phi(h*z)*(g3(z,1+r*model_v(z))-g3(z,1))/r,cutoff)
            weak.append({'n':n,'fractional_ratio':mp.nstr(result/expected,18)})
        q=lambda z:(1-4*z*z)**4
        bump_expected=phi(0)*mp.quad(lambda z: phi(z)*(z*z-1)*q(z)/2,[-mp.mpf('.5'),0,mp.mpf('.5')])
        pointwise=[]
        for n in (16,32,64,96):
            u=mp.power(2,-n)
            width=u**(alpha/2)
            amp=u**((alpha-2)/2)
            def bump(z):
                base=(1+width*z/u)**alpha
                return phi(u+width*z)*(p(z,base*(1+amp*q(z)))-p(z,base))/amp
            result=mp.quad(bump,[-mp.mpf('.5'),0,mp.mpf('.5')])
            pointwise.append({'n':n,'local_bump_ratio':mp.nstr(result/bump_expected,18)})
        return {'digits':digits,'alpha':str(alpha),'weak_coefficient_fixed_core':mp.nstr(expected,18),
                'weak_calibration_Gaussian_test':weak,'pointwise_bump_coefficient':mp.nstr(bump_expected,18),
                'pointwise_calibration_locally':pointwise}


print(json.dumps({'python':platform.python_version(),'sympy':sp.__version__,'mpmath':mp.__version__,
    'exact_checks':'three vanishing moments and integration by parts for alpha=5/2,7/2',
    'numerical_runs':[calibrate(digits) for digits in (50,75)]}, indent=2))
