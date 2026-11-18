"""Rationally certify an integer-difference witness for the R002 two-layer trace."""

import json
import sys
from fractions import Fraction as F
from math import factorial, isqrt
from pathlib import Path

import mpmath as mp


J,K=80,70
S=512
Q=10**27
DELTA=F(1,16)
CUT=F(1,2)
LAG=F(3,8)
W=F(2,3)


def floor(q):
    return q.numerator//q.denominator


def ceil(q):
    return -floor(-q)


def exp_neg(t):
    assert 0<=t<=1
    terms=[F(1)]
    for i in range(1,22):
        terms.append(-terms[-1]*t/i)
    upper=sum(terms[:-1])
    lower=upper+terms[-1]
    return lower,upper


def invsqrt_lower(t):
    return isqrt((Q*Q*t.denominator)//t.numerator)


def bounded_kernel(t):
    p=F(J,J+K)
    level=p**J*(1-p)**K
    elo,ehi=exp_neg(t)
    lo=floor(Q*(elo**J)*(1-ehi)**K/level)
    hi=min(Q,ceil(Q*(ehi**J)*(1-elo)**K/level))
    return lo,hi,elo,ehi


def rational_certificate():
    p=F(J,J+K)
    start=int(DELTA*S)
    assert F(start,S)==DELTA and S%2==0
    endpoints=[bounded_kernel(F(i,S)) for i in range(start,S+1)]
    assert endpoints[0][2]>p
    # K is increasing in t for t<=DELTA. Integrate t^-1/2 exactly there.
    assert 4*DELTA==F(1,4)
    sum_num=endpoints[0][1]*S*Q//2
    for i in range(start,S):
        lo,hi=F(i,S),F(i+1,S)
        e0,e1=endpoints[i-start],endpoints[i-start+1]
        lower_k=min(e0[0],e1[0])
        might_hit_peak=e1[2]<=p<=e0[3]
        upper_k=Q if might_hit_peak else max(e0[1],e1[1])
        q_upper=invsqrt_lower(lo)+1
        if lo>=CUT:
            q_upper-=(2*invsqrt_lower(hi-LAG))//3
        else:
            assert hi<=CUT
        sum_num+=(upper_k if q_upper>=0 else lower_k)*q_upper
    result=F(sum_num,S*Q*Q)
    assert result<F(-1,1000),result
    return result


def numerics(dps,method,j,k):
    with mp.workdps(dps):
        p=mp.mpf(j)/(j+k)
        log_peak=j*mp.log(p)+k*mp.log1p(-p)

        def kn(t):
            if t==0:
                return mp.mpf(0)
            return mp.exp(-j*t+k*mp.log(-mp.expm1(-t))-log_peak)

        free=2*mp.quad(lambda x: kn(x*x),(0,mp.mpf('.25'),mp.sqrt(mp.mpf('.5')),1),method=method)
        sub=mp.mpf(2)/3*mp.quad(lambda t: kn(t)/mp.sqrt(t-mp.mpf(3)/8),
                                (mp.mpf('.5'),mp.mpf(27)/40,1),method=method)
        return free-sub


def main():
    assert F(32,27)>F(13,12)**2
    xleft=(F(11,20))
    xright=(F(3,5))
    assert exp_neg(F(3,5))[1]<xleft
    assert exp_neg(F(1,2))[0]>xright
    assert 17600//1955>=9 and sum(F(9)**j/factorial(j) for j in range(15))>4416
    exact=rational_certificate()
    with mp.workdps(70):
        cal=[]
        for j,k in ((59,53),(J,K)):
            main=numerics(70,'tanh-sinh',j,k)
            check=numerics(50,'gauss-legendre',j,k)
            difference=abs(main-check)
            assert difference<mp.mpf('1e-44')
            assert main<0
            cal.append({'j':j,'k':k,'scaled_integral_without_gaussian_constant':mp.nstr(main,42),
                        'absolute_cross_method_difference':mp.nstr(difference,13)})
        cert={'scope':'Two-layer R002 model at y=0; rational finite-difference witness plus supplemental numerical calibration',
              'python':sys.version,'mpmath':mp.__version__,
              'rational_witness':{'j':J,'k':K,'grid':S,'fixed_point_denominator':str(Q),
                                  'normalized_integral_upper_bound_exact_fraction':str(exact),
                                  'normalized_integral_upper_bound_decimal':
                                      mp.nstr(mp.mpf(exact.numerator)/exact.denominator,40)},
              'generic_kl_witness':{'n':17600,'j':10120,'k':7480,
                                    'negative_t_interval_from_rational_x':[str(xleft),str(xright)],
                                    'lower_density_magnitude':'1/24','inner_negative_mass_lower':'1/2208'},
              'calibrations':cal}
    path=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).with_name('two-layer-results.json')
    path.write_text(json.dumps(cert,indent=2)+'\n')
    print(cert['rational_witness']['normalized_integral_upper_bound_decimal'])


if __name__=='__main__':
    main()
