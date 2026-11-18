"""Independent high-precision quadrature and exact two-layer moment signs."""
import json
import math
import platform
from fractions import Fraction as F
from pathlib import Path
import mpmath as m
import sympy as s

m.mp.dps = 70


def density(v,y):
    if not v:
        return m.mpf(0)
    return m.exp(-y*y/(2*v))/m.sqrt(2*m.pi*v)


def three_layer(j,k,y):
    half=m.mpf(1)/2
    e=3*m.exp(-m.mpf(16)/3)/5
    peak=(m.mpf(j)/(j+k))**j*(m.mpf(k)/(j+k))**k
    def w(t):
        return m.exp(-j*t)*(-m.expm1(-t))**k/peak
    ranges=([0,.25,.4,.48,.5],[.5,.52,.54,.57,.61,.7,1],[1,1.1,1.2,1.5])
    models=(lambda t:density(t,y),
            lambda t:density(t,y)-e*density(t-3*half/4,y-2),
            lambda t:e*density(t-7*half/4,y-2))
    segments=[m.quad(lambda t:w(t)*g(t), ts) for g,ts in zip(models,ranges)]
    return dict(j=j,k=k,y=m.nstr(y,50), normalized_segments=[m.nstr(x,55) for x in segments],
                normalized_density=m.nstr(sum(segments),55), raw_density=m.nstr(sum(segments)*peak,55))


def constant(y,q):
    a=m.sqrt(2*q)
    direct=m.quad(lambda t:m.exp(-q*t)*density(t,y),[0,.02,.1,.5,1])
    tail=lambda x:m.erfc(x/m.sqrt(2))/2
    closed=(m.exp(-a*y)*tail(y-a)-m.exp(a*y)*tail(y+a))/a
    assert abs(direct-closed)<m.mpf('1e-65')
    return [m.nstr(y,8),m.nstr(q,8),m.nstr(direct,55),m.nstr(abs(direct-closed),10)]


def second_model():
    def moment(n):
        g=F(math.prod(range(1,2*n,2)))
        return g*(F(1,2)**n-F(11,20)*F(1,8)**n)
    def M(k):
        return s.Matrix(k+1,k+1, lambda i,j:s.Rational(moment(i+j)))
    signs=[int(s.sign(M(k).det(method='bareiss'))) for k in range(14)]
    assert signs==[1]*13+[-1]
    h=M(13)
    witness=(-h[:13,:13].inv()*h[:13,13]).col_join(s.Matrix([1]))
    exact=(witness.T*h*witness)[0]
    assert exact<0
    return dict(leading_principal_minor_signs=signs,
                witness_even_coefficients=[str(x) for x in witness],
                exact_witness_expectation=str(exact),
                recovery_time='50/93',
                gaussian_marginal_variance='29/40',
                gaussian_marginal_fourth_moment='381/160',
                r2_tail_mass_w= '11/20')


def main():
    three=[three_layer(1000,709,m.mpf(y)) for y in ('2.829','2.8295','2.830')]
    assert all(m.mpf(x['normalized_density'])<0 for x in three)
    pairs=[('0.35','0.1'),('0.8','1.4'),('2','3')]
    result=dict(versions=[platform.python_version(),m.__version__,s.__version__],dps=m.mp.dps,
                three_layer=three, constant_control=[constant(m.mpf(y),m.mpf(q)) for y,q in pairs],
                two_layer=second_model())
    Path('independent.json').write_text(json.dumps(result,indent=2)+'\n')
    print(result['two_layer']['leading_principal_minor_signs'])
    for record in three:
        print(record['y'],record['normalized_density'],record['raw_density'])

if __name__=='__main__':
    main()
