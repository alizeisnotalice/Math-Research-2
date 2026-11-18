from pathlib import Path
import json

import sympy as s

ROOT=Path(__file__).parent


def main():
    k=s.symbols('kappa', positive=True)
    y=s.symbols('y', real=True)
    x=s.symbols('z', positive=True)
    # Gaussian termwise cosine expansion, normalized by its zeroth coefficient.
    ratios=[]
    for n in range(25):
        from_gauss=(-1)**n*2**n*s.rf(k+s.Rational(1,2),n)/s.factorial(2*n)
        from_hyper=(-s.Rational(1,2))**n*s.rf(k+s.Rational(1,2),n) / (
                       s.rf(s.Rational(1,2),n)*s.factorial(n))
        assert s.factor(from_gauss-from_hyper)==0
        ratios.append(str(s.factor(from_gauss)))
    v=s.Function('H')
    q=v(y**2/2)
    operator=s.diff(q,y,2)-y*s.diff(q,y)+2*k*q
    kom=2*(x*s.diff(v(x),x,2)+(s.Rational(1,2)-x)*s.diff(v(x),x)+k*v(x))
    assert s.simplify(operator-kom.subs(x,y**2/2))==0
    gauss_tail=[s.factor(2**n*s.rf(k+1,n)*s.rf(k+s.Rational(1,2),n)/s.factorial(n))
                for n in range(5)]
    expected=[s.rf(2*k+1,2*n)/(2**n*s.factorial(n)) for n in range(5)]
    assert all(s.factor(a-b)==0 for a,b in zip(gauss_tail,expected))
    # c_k*Gamma(1/2)/Gamma(-k), multiplied by the Marchaud factor
    # and the conversion from z^{-k-1/2} to |y|^{-2k-1}.
    c=2**(k-s.Rational(1,2))*s.gamma(k+s.Rational(1,2))/s.pi
    multiplier=s.gamma(1-k)/(k*2**k)
    tail=s.combsimp(c*s.sqrt(s.pi)*2**(k+s.Rational(1,2))*multiplier/s.gamma(-k))
    claim=-2**k*s.gamma(k+s.Rational(1,2))/s.sqrt(s.pi)
    assert s.simplify(tail-claim)==0
    # Small-kappa root after writing log(L/pi) as an independent formal symbol.
    u,w,b1,b2=s.symbols('u ell b1 b2')
    z=1/u+w/2+b1*u+b2*u**2
    residual=b1*u+b2*u**2-s.log(z*u)/2+1/(2*z)+s.Rational(5,8)/z**2
    expansion=s.series(residual,u,0,3).removeO().expand()
    solution=s.solve([expansion.coeff(u,1),expansion.coeff(u,2)],(b1,b2))
    assert s.simplify(solution[b1]-(w/4-s.Rational(1,2)))==0
    assert s.simplify(solution[b2]-(-w**2/16+3*w/8-s.Rational(7,8)))==0
    result=dict(sympy=s.__version__,cosine_series_ratios=ratios,
                marchaud_tail=str(tail),tail_expansion_coefficients=list(map(str,gauss_tail)),
                normalized_ode=str(kom),root_1_over_L=str(solution[b1]),root_1_over_L2=str(solution[b2]))
    (ROOT/'symbolic.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='cosine_series_ratios'}))


if __name__=='__main__':
    main()
