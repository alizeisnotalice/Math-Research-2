"""Finite exterior moment controls for a smoothed quartic core."""

import json
import sys
from pathlib import Path

import mpmath as mp

from jet_check import cutoff


TILTS = (0,-1,1,2)
R_LABEL = '0.0001'
BUMP_RADIUS = '0.12'


def run(dps, method):
    with mp.workdps(dps):
        R = mp.mpf(R_LABEL)
        radius = mp.mpf(BUMP_RADIUS)
        nodes = [mp.log(x) for x in (mp.mpf(1)/4,mp.mpf(1)/2,2,4)]
        q = 1/mp.sqrt(2*mp.pi)
        cuts = [-1,-mp.mpf('0.5'),0,mp.mpf('0.5'),1]

        def normal(x):
            return q*mp.exp(-x*x/2)

        def response(t,x,d):
            if t == 0:
                return normal(x)*d
            return normal(x)*mp.exp(t*x+t*t/2)*mp.expm1(t*t*d/2)

        core = [2*R*mp.quad(
            lambda u: response(t,2*R*u,cutoff(u,mp.mpf(1))*((2*R*u)**4-1)),
            cuts,method=method) for t in TILTS]

        norm = [radius*mp.quad(lambda u: normal(x+radius*u)*cutoff(u,mp.mpf(1)),
                               cuts,method=method) for x in nodes]

        def equations(*a):
            return [base+sum(radius*mp.quad(
                lambda u: response(t,x+radius*u,at*cutoff(u,mp.mpf(1))/z),
                cuts,method=method) for x,at,z in zip(nodes,a,norm)) for t,base in zip(TILTS,core)]

        def jacobian(*a):
            entries = []
            for t in TILTS:
                row = []
                for x,z,at in zip(nodes,norm,a):
                    if t == 0:
                        row.append(mp.mpf(1))
                    else:
                        def diff(u):
                            y = x+radius*u
                            psi = cutoff(u,mp.mpf(1))/z
                            return normal(y)*mp.exp(t*y+t*t*(1+at*psi)/2)*t*t*psi/2
                        row.append(radius*mp.quad(diff,cuts,method=method))
                entries.append(row)
            return mp.matrix(entries)

        j0 = jacobian(0,0,0,0)
        start = mp.lu_solve(j0, -mp.matrix(core))
        answer = mp.findroot(equations,tuple(start),J=jacobian)
        residuals = equations(*answer)
        height = max(abs(a/z) for a,z in zip(answer,norm))
        assert height<mp.mpf('0.5')
        assert max(map(abs,residuals))<mp.mpf(10)**(-dps+5)
        return mp.det(j0), tuple(answer), tuple(residuals), height, tuple(core)


def main():
    det, params, residuals, height, core = run(54,'gauss-legendre')
    _, check, *_ = run(38,'tanh-sinh')
    with mp.workdps(54):
        error = max(abs(a-b) for a,b in zip(params,check))
        assert error<mp.mpf('1e-30')
        result = {'scope':'Numerical exterior inverse-function check for a quartic core, with no all-jet numerical claim.',
                  'python':sys.version,'mpmath':mp.__version__,
                  'tilts_including_0_for_variance':TILTS, 'b':1,'R':R_LABEL,
                  'nodes':'log(1/4),log(1/2),log(2),log(4)',
                  'bump_radius':BUMP_RADIUS, 'jacobian_determinant':mp.nstr(det,35),
                  'core_residual':[mp.nstr(a,30) for a in core],
                  'parameter_solution':[mp.nstr(a,35) for a in params],
                  'matched_residual':[mp.nstr(a,15) for a in residuals],
                  'max_exterior_variance_deviation_upper_bound':mp.nstr(height,30),
                  'cross_method_parameter_error':mp.nstr(error,12),
                  'methods':{'gauss-legendre':54,'tanh-sinh':38}}
    path = Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).with_name('exterior-results.json')
    path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['jacobian_determinant'],result['cross_method_parameter_error'])


if __name__=='__main__':
    main()
