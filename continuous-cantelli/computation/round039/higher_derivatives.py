"""Calibrate absolute translate derivative bounds for a one-sided D envelope."""

import itertools
import json
import sys
from pathlib import Path

import mpmath as mp
import sympy as sp


ORDERS = (2, 3, 4, 6, 8)
CASES = ('plain:1000', 'plain:1000000', 'cusp:500', 'cusp:1000000')


def coefficients(n):
    return tuple(sp.Rational((-1)**j*sp.factorial(n), 2**j*sp.factorial(j)*sp.factorial(n-2*j))
                 for j in range(n//2+1))


def exact_checks():
    W, f, s, x = sp.symbols('W f s x')
    output = {}
    for n in ORDERS:
        coeffs = coefficients(n)
        expected = sum(c*f**j*W**(n-2*j) for j,c in enumerate(coeffs))
        generated = sp.series(sp.exp(s*W-s*s*f/2), s, 0, n+1).removeO().coeff(s,n)*sp.factorial(n)
        assert sp.expand(generated-expected) == 0
        # E[e^{sG} He_n(G)] / E[e^{sG}] = E[He_n(G+s)].
        raw_moment = lambda k: 0 if k%2 else sp.factorial2(k-1)
        shifted_hermite = sp.Poly(sp.hermite_prob(n,x+s),x)
        transform = sum(c*raw_moment(k) for (k,),c in shifted_hermite.terms())
        assert sp.expand(transform-s**n) == 0
        output[n] = [str(c) for c in coeffs]
    return output


def quadrature(label, variant, dps, method):
    with mp.workdps(dps):
        name, count = label.split(':')
        if name == 'cusp':
            center = mp.pi*int(count)
            L = mp.findroot(lambda x: x-mp.log(2*x)-center, center+mp.log(2*center))
        else:
            L = mp.mpf(count)
        root = mp.sqrt(L)
        base = L-mp.log(2*L)
        common = mp.log(2)/2-mp.mpf('0.1')
        downward_at_center = mp.mpf('0.1') if variant == 'bounded' else mp.sqrt(base)
        penalty_at_center = downward_at_center*abs(mp.sin(base))
        mass_scale = mp.exp(-penalty_at_center)
        low = mp.exp(-2*L)
        cuts = [-15, -9, -5, 0, 5, 9, 15]
        for index in range(int(mp.floor(base/mp.pi))-1, int(mp.ceil(base/mp.pi))+2):
            cusp = root*mp.expm1(base-index*mp.pi)
            if -15<cusp<15 and min(abs(cusp-c) for c in cuts)>mp.mpf('1e-30'):
                cuts.append(cusp)
        cuts.sort()

        def centered(z):
            ell = base-mp.log1p(z/root)
            y = 2*L+2*root*z
            downward = mp.mpf('0.1') if variant == 'bounded' else mp.sqrt(ell)
            delta = -mp.log(y)+mp.log(ell)/2+common-downward*abs(mp.sin(ell))
            D = L+delta
            gap = 2*root*z-2*delta
            w = -y*gap/(2*D*root)
            f = y*y/(2*D*L)
            density = 2*root/mp.sqrt(2*mp.pi)*mp.exp(delta+penalty_at_center-gap*gap/(4*D)-low*y*y/2)
            return density, w, f

        saddle = mp.findroot(lambda z: centered(z)[1], (-1, 1))
        if min(abs(saddle-c) for c in cuts)>mp.mpf('1e-30'):
            cuts.append(saddle)
            cuts.sort()

        def absolute(z, n, coeffs):
            density, w, f = centered(z)
            return density*sum(c*f**j*abs(w)**(n-2*j) for j,c in enumerate(coeffs))

        integrals = {}
        for n in ORDERS:
            coeffs = tuple(abs(mp.mpf(int(c.p))/int(c.q)) for c in coefficients(n))
            integrals[n] = mass_scale*mp.quad(lambda z: absolute(z,n,coeffs), cuts, method=method)
        mass = mass_scale*mp.quad(lambda z: centered(z)[0], cuts, method=method)
        return L, mass, integrals


def main():
    exact = exact_checks()
    rows = []
    with mp.workdps(72):
        for label, variant in itertools.product(CASES, ('bounded','downward_sqrt')):
            L, mass, primary = quadrature(label, variant, 72, 'gauss-legendre')
            _, alt_mass, check = quadrature(label, variant, 48, 'tanh-sinh')
            difference = max(abs(alt_mass-mass), *(abs(primary[n]-check[n]) for n in ORDERS))
            relative = max(abs(alt_mass-mass)/mass, *(abs(primary[n]-check[n])/primary[n] for n in ORDERS))
            assert difference<mp.mpf('1e-30'), difference
            assert relative<mp.mpf('1e-30'), relative
            rows.append({'case':label,'downward_penalty':variant,'L':mp.nstr(L,30),
                         'mass':mp.nstr(mass,32),
                         'absolute_derivative_term_sum_times_t_over_sqrt_L_power_n':
                             {n:mp.nstr(v,30) for n,v in primary.items()},
                         'cross_method_max_difference':mp.nstr(difference,12),
                         'cross_method_max_relative_difference':mp.nstr(relative,12)})
            print(label, variant, mp.nstr(mass,10), mp.nstr(primary[2],10), mp.nstr(primary[8],10), flush=True)
    output = {'scope':'Numerical centered origin model integrals and symbolic derivative/Hermite check',
              'python':sys.version,'mpmath':mp.__version__,'sympy':sp.__version__,
              'methods_digits':{'gauss-legendre':72,'tanh-sinh':48},
              'kernel_coefficients_exact':exact,'rows':rows}
    path = Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).with_name('higher-derivative-results.json')
    path.write_text(json.dumps(output,indent=2)+'\n')


if __name__=='__main__':
    main()
