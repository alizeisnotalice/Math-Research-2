"""Exact triangular algebra and a numerical two-equation local potential model."""

import json
import sys
from pathlib import Path

import mpmath as mp
import sympy as sp


u, c, a, b = sp.symbols('u c a b')
jet = sp.symbols('w1:7')
g_relative = sp.series(sp.exp(-u-u*u/2), u, 0, 7).removeO()


def diagonal(n, w, g):
    return sp.expand(sum(
        sp.factorial(n+2*j-2)/(2**j*sp.factorial(j)*sp.factorial(n))
        * (g*w**j).expand().coeff(u,n-2*j-2)
        for j in range(1,(n-2)//2+1)
    ))


def verify_triangle():
    w = c+sum(value*u**j for j,value in enumerate(jet,1))
    rows = {}
    for n in range(5,11):
        value = diagonal(n,w,g_relative)
        k = n-4
        assert sp.diff(value, jet[k-1]) == sp.Rational(1,2)
        assert not any(value.has(x) for x in jet[k:])
        rows[n] = str(value)
    return rows


def cutoff(u, radius):
    q = abs(u)/radius
    if q <= mp.mpf('0.5'):
        return mp.mpf(1)
    if q >= 1:
        return mp.mpf(0)
    left = mp.exp(-1/(q-mp.mpf('0.5')))
    right = mp.exp(-1/(1-q))
    return right/(left+right)


def finite_example(dps, method):
    with mp.workdps(dps):
        n_values = (6,8)
        w_formal = 1+a*u*u+b*u**4
        source_formal = 1-u*u/2+u**4/8
        polys = tuple(diagonal(n,w_formal,source_formal) for n in n_values)
        expressions = tuple(sp.lambdify((a,b),p,'mpmath') for p in polys)
        r = mp.mpf('0.12')
        eps_a = mp.mpf('0.04')
        eps_b = mp.mpf('0.025')
        cuts = [mp.mpf(x) for x in ('0','0.0125','0.02','0.025','0.04','0.06','0.09','0.12')]
        hermite = {
            n: tuple(map(int,sp.Poly(sp.hermite_prob(n-2,u),u).all_coeffs()))
            for n in n_values
        }
        q0 = 1/mp.sqrt(2*mp.pi)

        def equations(a2,b4):
            values = []
            for n, expr in zip(n_values,expressions):
                def integrand(x):
                    w = 1+a2*x*x*cutoff(x,eps_a)+b4*x**4*cutoff(x,eps_b)
                    s = x**4*w
                    z = x/mp.sqrt(s)
                    return q0*mp.polyval(hermite[n],z)*mp.exp(-(z*z+x*x)/2)/s**(mp.mpf(n-1)/2)
                integral = 2*mp.quad(integrand,cuts,method=method)/mp.factorial(n)
                values.append(expr(a2,b4)+integral)
            return values

        formal_solution = sp.solve(polys,(a,b))
        if isinstance(formal_solution,list):
            formal_solution = formal_solution[0]
        else:
            formal_solution = tuple(formal_solution[q] for q in (a,b))
        start = tuple(mp.mpf(int(x.p))/int(x.q) for x in formal_solution)
        sol = mp.findroot(equations,start)
        residual = equations(*sol)
        bound = abs(sol[0])*eps_a**2+abs(sol[1])*eps_b**4
        assert bound < mp.mpf('0.5')
        assert max(map(abs,residual)) < mp.mpf(10)**(-dps+6)
        return tuple(sol), tuple(residual), bound, tuple(map(str,polys)), tuple(map(str,formal_solution)), str(r)


def main():
    matrix = verify_triangle()
    primary = finite_example(60,'gauss-legendre')
    check = finite_example(44,'tanh-sinh')
    with mp.workdps(60):
        error = max(abs(x-y) for x,y in zip(primary[0],check[0]))
        assert error < mp.mpf('1e-35'), error
        out = {'scope':'Symbolic audit of local triangular coefficients and a finite numerical model; no full-jet numerical claim.',
               'python':sys.version,'sympy':sp.__version__,'mpmath':mp.__version__,
               'symbolic_diagonal_by_order_relative_to_source_at_zero':matrix,
               'zero_linear_channel_by_even_symmetry':[5,7],
               'solved_orders':[6,8],
               'methods':{'gauss-legendre':60,'tanh-sinh':44},
               'cutoff_radius':primary[5], 'borel_radii':['0.04','0.025'],
               'diagonal_polynomials':primary[3],'formal_without_integral':primary[4],
               'numeric_a2_a4':[mp.nstr(v,50) for v in primary[0]],
               'non_diagonal_correction':[mp.nstr(v-mp.mpf(str(p)),30) if '/' not in p else mp.nstr(v-mp.mpf(p),30)
                                          for v,p in zip(primary[0],primary[4])],
               'actual_coefficient_residual_relative_to_g_at_zero':[mp.nstr(v,12) for v in primary[1]],
               'sup_norm_departure_from_one_upper_bound':mp.nstr(primary[2],30),
               'cross_method_max_difference':mp.nstr(error,12)}
    path = Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).with_name('jet-results.json')
    path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['numeric_a2_a4'], out['cross_method_max_difference'])


if __name__=='__main__':
    main()
