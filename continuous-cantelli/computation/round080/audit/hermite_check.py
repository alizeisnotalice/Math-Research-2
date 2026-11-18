import json
import math
import platform

import mpmath as mp
import sympy as sp

x,z,b=sp.symbols('x z b')
N=14
raw=sp.Poly(
    sp.expand(sum(x**i*z**i/sp.factorial(i) for i in range(N+1))
              *sum((b-1)**k*z**(2*k)/(2**k*sp.factorial(k)) for k in range(N//2+1))),
    z,b,
)
amount=0
for n in range(N-1):
    for j in range(1,(n+2)//2+1):
        m=n+2-2*j
        expected=sp.hermite_prob(m,x)/(2**j*sp.factorial(j)*sp.factorial(m))
        assert sp.expand(raw.coeff_monomial(z**(n+2)*b**j)-expected)==0
        amount+=1


def calibrate(digits):
    with mp.workdps(digits):
        def G(u):
            return (mp.expm1(u)-u)/u if abs(u)>mp.mpf('1e-10') else mp.nsum(lambda k:u**k/mp.factorial(k+1),[1,mp.inf])
        u=mp.findroot(lambda u:mp.expm1(u)-2*u,(1,2))
        c=2*u
        exps=[]
        for s in map(mp.mpf,('0.25','0.5','1','2')):
            moments=c**s*mp.quad(lambda q:1-G(c*q**(1/s)/2)**2,[0,mp.mpf('.5'),1])
            exps.append({'s':str(s),'layer_cake_lower_bound':mp.nstr(moments,18)})
        def laplace(time):
            return 1-2*time*mp.quad(lambda x:mp.exp(-2*time*x)*(1-G(x/2)**2),[0,1,c])
        return {'digits':digits,'c_star':mp.nstr(c,18),'fractional_moments':exps,
                'OU_time_epsilon_over_2_L2_squared_bound':mp.nstr(laplace(mp.mpf('.5')),18),
                'direct_OU_norm_contraction':mp.nstr(mp.e-2,18),
                'direct_OU_derivative_constant':mp.nstr(-2*mp.log(mp.e-2),18)}

partials=[]
partial=0.0
for k in range(1,8193):
    partial+=math.exp(k*(math.log(k)-1)-math.lgamma(k+2))
    if k in (64,256,1024,8192):
        partials.append({'count':k,'binary64_partial_tree_sum':partial})

print(json.dumps({'python':platform.python_version(),'sympy':sp.__version__,'mpmath':mp.__version__,
    'exact_symbolic_count':amount,'max_raw_transform_degree':N,
    'mpmath_calibration':[calibrate(50),calibrate(75)],'tree_sum_partial_values':partials},indent=2))
