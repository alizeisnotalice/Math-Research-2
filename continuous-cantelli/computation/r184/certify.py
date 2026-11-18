"""Exact Fraction certification of zG-G' crossings by positive moment recurrences."""
import argparse,json,os,platform,time
from fractions import Fraction as Q
from pathlib import Path
import mpmath as m


M=512
GAMMA_ORDER=2048
SERIES_ORDER=208
DIGITS=90


def gamma_normalized(a):
    """Enclose Gamma(a)/M**a by rational bounds for 0<a<1."""
    assert 0<a<1 and GAMMA_ORDER%2==0
    power=Q(1);s=power/a
    for n in range(1,GAMMA_ORDER+1):
        power*=Q(-M,n)
        term=power/(a+n)
        s+=term
    # Taylor integrated on [0,M]; the normalized positive tail is <= e^-M/M < 2^-M/M.
    return s-term,s+Q(1,M*2**M)


def base_ratio(k):
    a=gamma_normalized((1-k)/2);b=gamma_normalized((2-k)/2)
    assert a[0]>0 and b[0]>0 and 2*M==32**2
    return 32*b[0]/a[1],32*b[1]/a[0]


def parity_sums(k,z):
    """Enclose the even and odd sums; the odd sum omits the base multiplier mu1/mu0."""
    sums=[]
    for parity,term in [(0,z),(1,z*z/2)]:
        s=0*z
        for n in range(parity,SERIES_ORDER,2):
            s+=term
            term*=((n+1-k)*z*z)/((n+2)*(n+3))
        rho=z*z/(SERIES_ORDER+3)
        assert rho<1
        sums.append((s,s+term/(1-rho)))
    return sums


def interval(k,z,rat):
    even,odd=parity_sums(k,z)
    rlo,rhi=rat
    return -rhi+k*(even[0]+rlo*odd[0]),-rlo+k*(even[1]+rhi*odd[1])


def fmt_grid(n,d):
    s='-' if n<0 else ''
    x,y=divmod(abs(n),10**d)
    return f'{s}{x}.{y:0{d}d}'


def outward(pair,d):
    lo,hi=pair;s=10**d
    return [fmt_grid(lo.numerator*s//lo.denominator,d),fmt_grid(-((-hi.numerator*s)//hi.denominator),d)]


def mpq(r):return m.mpf(r.numerator)/r.denominator


def solve(k):
    ratio=base_ratio(k)
    with m.workdps(182):
        kk=mpq(k);a0=(1-kk)/2;a1=(2-kk)/2
        r=m.sqrt(2)*m.gamma(a1)/m.gamma(a0)
        assert mpq(ratio[0])<=r<=mpq(ratio[1])
    with m.workdps(122):
        kk=mpq(k);a0=(1-kk)/2
        r=+r
        def approx(z):
            e,o=parity_sums(kk,z)
            return -r+kk*(e[0]+r*o[0])
        root=m.findroot(approx,(m.mpf(0),r/kk),tol=m.mpf('1e-118'))
        assert 0<root<2
        den=10**DIGITS
        j=int(m.floor(root*den))
        u,v=Q(j-1,den),Q(j+2,den)
        left,right=interval(k,u,ratio),interval(k,v,ratio)
        assert 0<u<v and left[1]<0<right[0]
        def direct_integrand(y):
            t=y**(1/(1-kk))
            return (root-t)*m.exp(root*t-t*t/2)/(1-kk)
        mesh=[m.mpf(t)**(1-kk) for t in (0,.03125,.125,.5,1,2,3,5,8,12,20,40)]
        d_direct=m.quad(direct_integrand,mesh,method='tanh-sinh',maxdegree=9)/(2**(a0-1)*m.gamma(a0))
        mp_residual=abs(d_direct-approx(root))
        assert mp_residual<m.mpf('1e-105')
        max_sampling_tail=max(mpq(y[1]-y[0]) for z in (u,v) for y in parity_sums(k,z))
        return dict(kappa=str(k),certified_root_open=[fmt_grid(j-1,DIGITS),fmt_grid(j+2,DIGITS)],
                    normalized_D_at_left=outward(left,DIGITS+12),normalized_D_at_right=outward(right,DIGITS+12),
                    initial_moment_ratio=outward(ratio,DIGITS+12),approximate_root=m.nstr(root,110),
                    normalized_direct_comparison=m.nstr(mp_residual,35),
                    exact_sum_tail_upper_display=m.nstr(max_sampling_tail,35),
                    base_ratio_width_display=m.nstr(mpq(ratio[1]-ratio[0]),35))


def main():
    p=argparse.ArgumentParser();p.add_argument('kappa');p.add_argument('out',type=Path);args=p.parse_args()
    t=time.time();cpu=time.process_time();res=solve(Q(args.kappa))
    res.update(pid=os.getpid(),start_epoch=t,end_epoch=time.time(),cpu_seconds=time.process_time()-cpu,
        python=platform.python_version(),mpmath=m.__version__,mp_dps=122,gamma_order=GAMMA_ORDER,series_order=SERIES_ORDER)
    args.out.write_text(json.dumps(res,indent=2)+'\n')

if __name__=='__main__':main()
