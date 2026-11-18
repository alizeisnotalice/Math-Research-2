"""Precision stress tests of the early-interval heat cost and its upper-time inverse."""
import argparse,json,os,platform,time
from pathlib import Path
import mpmath as m


def cut_quad(fun,limit):
    nodes=[m.mpf(x) for x in (0,.25,1,3,7,15,31,63,127,256) if x<limit]
    return m.quad(fun,nodes+[limit],method='tanh-sinh',maxdegree=8)


def series_h(r):
    rr=r*r
    term=m.mpf(1);total=term
    for n in range(1,45):
        term*=-(2*n+1)/rr
        total+=term
    tail=abs(term)*(2*45+1)/rr
    return total,tail


def erfc_h(r):
    return r*r*(1-r*m.sqrt(m.pi/2)*m.exp(r*r/2)*m.erfc(r/m.sqrt(2)))


def normalized_heat(x,b,width,h):
    def integrand(w):
        r=x-b+w/x
        return m.exp(-(1-2*b/x)*w-(w/x)**2)*(x/r)**2*h(r)
    limit=min(256,x*width)
    return cut_quad(integrand,limit)


def invert(x,log_a,b):
    lam=x*x/2
    target=log_a+m.log(lam)
    if b>0:
        def log_scaled(R):
            return m.log(cut_quad(lambda w:m.exp(-w)*(1-(R-w)/lam)**(-m.mpf('1.5')),min(256,R)))
        R=m.findroot(lambda R:R+log_scaled(R)-target,(target-m.mpf('.1'),target+m.mpf('.1')),tol=m.mpf('1e-84'))
    else:
        factor=m.exp(target)
        def equation(k):
            R=factor*k
            return k*m.quad(lambda u:m.exp(R*u)*(1-R*u/lam)**(-m.mpf('1.5')),[0,.5,1])-1
        R=factor*m.findroot(equation,(m.mpf('.99'),m.mpf('1.01')),tol=m.mpf('1e-90'))
    d=R/(lam-R)
    if b>0:
        def integrand(w):
            t=1+d*(1-w/R)
            return m.exp(-w/t)/m.sqrt(t)
        check=m.log(d/R)+R+m.log(cut_quad(integrand,min(256,R)))
    else:
        check=m.log(d)+m.log(m.quad(lambda u:m.exp(lam*d*u/(1+d*u))/m.sqrt(1+d*u),[0,.5,1]))
    assert abs(check-log_a)<m.mpf('1e-70')
    return d,R,check-log_a


def calculate(left,right,wing):
    with m.workdps(100):
        a,b=map(m.mpf,(left,right))
        if wing<0:a,b=-b,-a
        rows=[]
        constant=m.sqrt(2/m.pi)*m.exp(-b*b)
        for x in map(m.mpf,(128,4096,1048576)):
            main=normalized_heat(x,b,b-a,lambda r:series_h(r)[0])
            independent=normalized_heat(x,b,b-a,erfc_h)
            err=abs(main-independent)
            assert err<m.mpf('1e-68')
            log_a=m.log(constant)+b*x-3*m.log(x)+m.log(main)
            delay,R,residual=invert(x,log_a,b)
            g_series=1+4*b/x+(11*b*b-7)/x**2
            if b>0:
                time_series=2*b/x+(-2*m.log(x)+2*b*b-m.log(2*m.pi))/x**2
            elif b==0:
                time_series=constant/x**3-constant**2/(4*x**4)
            else:
                time_series=constant*m.exp(b*x)/x**3*g_series
            row=dict(X=int(x),beta=b,normalized_B=main,normalized_B_check=independent,
                     method_difference=err,log_B_over_g=log_a,B_over_g=m.exp(log_a),
                     endpoint_series=g_series,B_series_error=main-g_series,delay=delay,
                     delay_over_B_over_g=delay*m.exp(-log_a),normalized_R=R,
                     time_series=time_series,time_relative_error=delay/time_series-1,
                     inversion_log_residual=residual)
            rows.append({k:m.nstr(v,66) if isinstance(v,m.mpf) else v for k,v in row.items()})
    return rows


def main():
    p=argparse.ArgumentParser();p.add_argument('a');p.add_argument('b');p.add_argument('wing',type=int,choices=(-1,1));p.add_argument('out',type=Path);a=p.parse_args()
    start=time.time();cpu=time.process_time()
    rows=calculate(a.a,a.b,a.wing)
    info=dict(interval=[a.a,a.b],wing=a.wing,rows=rows,pid=os.getpid(),python=platform.python_version(),mpmath=m.__version__,
              start_epoch=start,end_epoch=time.time(),cpu_seconds=time.process_time()-cpu,dps=100)
    a.out.write_text(json.dumps(info,indent=2)+'\n')

if __name__=='__main__':main()
