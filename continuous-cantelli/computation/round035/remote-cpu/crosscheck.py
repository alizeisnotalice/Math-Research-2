"""Independent arbitrary-precision direct and Gaussian-dual quadrature checks."""
import json
import platform
from pathlib import Path

import mpmath as mp


def phi(x):
    return mp.exp(-x*x/2)/mp.sqrt(2*mp.pi)


def cdf(x):
    return mp.erfc(-x/mp.sqrt(2))/2


def interval_cdf(a,b):
    if a>0:
        return cdf(-a)-cdf(-b)
    return cdf(b)-cdf(a)


def gaussian_dual(lo,hi,points,f):
    def integrand(x):
        variance=f(x)
        root=mp.sqrt(variance)
        d1,d2=lo*root-x/root,hi*root-x/root
        action=x*x*(1-variance)/(2*variance)
        return mp.exp(action)*interval_cdf(d1,d2)/root
    return mp.quad(integrand,points)/(hi-lo)


def gpu_check(a,w,p,h,offset,half):
    a,w,p,h,offset,half=map(mp.mpf,[a,w,p,h,offset,half])
    lo,hi=a+offset-half,a+offset+half
    y_pts=[-w,-w/2,mp.mpf(0),w/2,w]
    fn=lambda x:1-h/a**2*(1-((x-a)/w)**2)**p
    primal=lambda s:mp.quad(
        lambda y:phi(y-s)*mp.exp(h*(1-(y/w)**2)**p*(1+s/a)**2/2),y_pts)
    p_integral=mp.quad(primal,[lo-a,offset,hi-a])/(hi-lo)
    d_integral=gaussian_dual(lo,hi,[a+y for y in y_pts],fn)
    assert abs(p_integral/d_integral-1)<mp.mpf('1e-40')
    return dict(a=str(a),w=str(w),p=str(p),h=str(h),offset=str(offset),halfwidth=str(half),
                primal=str(p_integral),dual=str(d_integral),relative_disagreement=str(abs(p_integral/d_integral-1)))


def mixed_check():
    t=mp.mpf(4096)
    lo,hi=t+mp.mpf('3.5'),t+mp.mpf('4.5')
    def D(x):
        ell=-mp.log(x)
        return ell+mp.log(2*ell)/2
    core_points=[mp.mpf(y)/t for y in (10,13,16)]
    origin=lambda x:x*x/(2*D(x))
    def core(s):
        return mp.quad(lambda x:phi(x)*mp.exp(s*x-s*s*origin(x)/2),core_points)
    def tail(s):
        return interval_cdf(t+2-s,t+6-s)
    p_integral=mp.quad(lambda s:core(s)+tail(s),[lo,t+4,hi])/(hi-lo)
    d_core=gaussian_dual(lo,hi,core_points,origin)
    d_tail=gaussian_dual(lo,hi,[t+2,t+4,t+6],lambda x:mp.mpf(1))
    assert abs(p_integral-d_core-d_tail)<mp.mpf('1e-40')
    return dict(origin_dual=str(d_core),tail_dual=str(d_tail),combined_dual=str(d_core+d_tail),
                combined_primal=str(p_integral),combined_at_nominal_t=str(core(t)+tail(t)))


def main():
    mp.mp.dps=58
    out=dict(python=platform.python_version(),mpmath=mp.__version__,digits=mp.mp.dps,
             gpu_top=gpu_check(8,6,1,12,4,mp.mpf(1)/32),mixed_finite_certificate=mixed_check())
    Path('crosscheck.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    main()
