"""Exact rational Gaussian-bridge endpoint enclosures for the fixed fat Cantor set."""
import argparse,json,math,os,platform,time
from fractions import Fraction as Q
from pathlib import Path
import mpmath as mp


S=Q(7,8);C=Q(15,16);LENGTH=Q(3,10000)
LEFT=Q(40001,20000);RIGHT=Q(40007,20000)
LO=Q(1139039,10**7);HI=Q(1166260,10**7)
GAP=Q(1,256)


def arctan(x,n):
    total=sum((-1)**i*x**(2*i+1)/(2*i+1) for i in range(n))
    return total,total+x**(2*n+1)/(2*n+1)


def pi_bounds(n):
    a,b=arctan(Q(1,5),n),arctan(Q(1,239),n)
    return 16*a[0]-4*b[1],16*a[1]-4*b[0]


def sqrt_bounds(x):
    scale=10**54
    k=math.isqrt((x.numerator*scale**2)//x.denominator)
    lo,hi=Q(k,scale),Q(k+1,scale)
    assert lo**2<=x<hi**2
    return lo,hi


def exp_bounds(x,n):
    s=sum(x**i/math.factorial(i) for i in range(n+1))
    return s,s+x**(n+1)/math.factorial(n+1)/(1-x/(n+2))


def interval_numerator(theta,n):
    # Integral from -1/2 to 1/2 of exp(-x*x/(2*theta)).
    # Even partial sum is an upper bound and the previous odd sum is a lower bound.
    assert n%2==0
    terms=[(-1)**i/(Q(8)*theta)**i/((2*i+1)*math.factorial(i)) for i in range(n+1)]
    return sum(terms[:-1]),sum(terms)


def w_bound(x,n,c):
    a,b=sqrt_bounds(S);p,q=exp_bounds(x*x/14,n)
    return c*a*p,c*b*q


def convert(x):return mp.mpf(x.numerator)/x.denominator


def verify(n,c,lo,hi):
    lower_w=w_bound(LEFT,n,c);upper_w=w_bound(RIGHT,n,c)
    pi=pi_bounds(n)
    assert n>=16 and n%2==0 and lower_w[0]>1
    rows=[]
    for xi in (lo,hi):
        theta=S*xi/(S+xi*LENGTH**2)
        a,b=interval_numerator(theta,n)
        sqrt2pi_lo=sqrt_bounds(2*pi[0]*theta)[0]
        sqrt2pi_hi=sqrt_bounds(2*pi[1]*theta)[1]
        assert a>GAP
        lower=(lower_w[0]*(a-GAP)/sqrt2pi_hi,lower_w[1]*(b-GAP)/sqrt2pi_lo)
        upper=(upper_w[0]*a/sqrt2pi_hi,upper_w[1]*b/sqrt2pi_lo)
        if xi==lo: assert lower[0]>1
        if xi==hi: assert upper[1]<1
        rxi=convert(xi);th=convert(theta)
        erf=mp.erf(1/mp.sqrt(8*th))
        alternate=mp.quad(lambda x:mp.exp(-x*x/(2*th)),[-mp.mpf('.5'),0,mp.mpf('.5')])/mp.sqrt(2*mp.pi*th)
        assert abs(erf-alternate)<mp.mpf('1e-67')
        f=convert(c)*mp.sqrt(convert(S))*mp.exp(convert(LEFT)**2/14)*(erf-convert(GAP)/mp.sqrt(2*mp.pi*th))
        g=convert(c)*mp.sqrt(convert(S))*mp.exp(convert(RIGHT)**2/14)*erf
        assert convert(lower[0])<=f<=convert(lower[1]) and convert(upper[0])<=g<=convert(upper[1])
        rows.append(dict(xi=str(xi),theta=str(theta),u=str(xi*LENGTH**2),t=str(S+xi*LENGTH**2),
                         geometric_lower_exact=list(map(str,lower)),geometric_upper_exact=list(map(str,upper)),
                         lower_ratio_display=mp.nstr(f,45),upper_ratio_display=mp.nstr(g,45),
                         erf_direct_quad_difference=mp.nstr(abs(erf-alternate),25)))
    assert hi*LENGTH**2<1-S
    return rows


def main():
    p=argparse.ArgumentParser();p.add_argument('terms',type=int);p.add_argument('out',type=Path);a=p.parse_args()
    t=time.time();cpu=time.process_time()
    with mp.workdps(78):
        rows=dict(mixed=verify(a.terms,C,LO,HI),deterministic=verify(a.terms,Q(1),Q(1472570,10**7),Q(1500079,10**7)))
    rec=dict(order=a.terms,rows=rows,pid=os.getpid(),python=platform.python_version(),mpmath=mp.__version__,
             start=t,end=time.time(),cpu_seconds=time.process_time()-cpu)
    a.out.write_text(json.dumps(rec,indent=2)+'\n')

if __name__=='__main__':main()
