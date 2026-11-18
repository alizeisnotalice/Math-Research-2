"""Independent Gaussian-time Marchaud quadratures and asymptotic coefficients."""
import argparse,json,os,platform,time
from fractions import Fraction as Q
from pathlib import Path
import mpmath as m


def mpq(x):return m.mpf(x.numerator)/x.denominator


def coefficients(k,v,n):
    values=[Q(1)]
    for j in range(1,n+1):
        values.append(values[-1]*(k+j)*(k+j-Q(1,2))*2*v/j)
    return values


def evaluate(kr,vr):
    with m.workdps(100):
        k,v=mpq(kr),mpq(vr)
        a=k+m.mpf('.5')
        g=m.gamma(a)
        cs=coefficients(kr,vr,4)
        c=[mpq(z) for z in cs]
        rows=[]
        for X in (16,64,1024,65536):
            x=m.mpf(X);xx=x*x
            p0=m.exp(-xx/(2*v))/m.sqrt(2*m.pi*v)
            norm=2**k/(m.sqrt(m.pi)*x**(1+2*k))
            def geth(y):return v*y**(1/(1-k))
            def direct_near(y):
                h=geth(y)
                if h==0:return -p0*(xx-v)/(2*v*v)
                exponent=-m.log1p(h/v)/2+xx*h/(2*v*(v+h))
                return -p0*m.expm1(exponent)/h
            def derivative_near(y):
                t=v+geth(y)
                return m.exp(-xx/(2*t))*(xx-t)/(2*t*t*m.sqrt(2*m.pi*t))
            jac=v**(1-k)/(1-k)
            near_direct=jac*m.quad(direct_near,[0,.0625,.25,.5,.75,1],method='tanh-sinh',maxdegree=9)
            near_parts=jac*m.quad(derivative_near,[0,.0625,.25,.5,.75,1],method='tanh-sinh',maxdegree=9)
            U=min(m.mpf(256),xx/(4*v));ar=2*v/xx
            mesh=[m.mpf(z)**a for z in (0,.03125,.125,.5,1,2,4,8,16,32,64,128) if z<U]+[U**a]
            def trans(y):return y**(1/a)
            def far_direct(y):
                r=trans(y)
                return m.exp(-r)*(1-ar*r)**(-1-k)/a
            def far_parts(y):
                r=trans(y)
                return m.exp(-r)*(2*r-1)*(1-ar*r)**(-k)/a
            # The two formulae differ before quadrature; in the second a derivative of p_t is integrated.
            d=m.quad(far_direct,mesh,method='tanh-sinh',maxdegree=9)
            # At r=0 use a logarithmic coordinate so Gaussian quadrature sees a smooth integrand.
            r0=m.mpf('.125')
            def parts_log(y):
                r=m.exp(-y)
                return m.exp(-r-a*y)*(2*r-1)*(1-ar*r)**(-k)
            log_mesh=[-m.log(r0),4,8,16,32,64,128,256,512]
            tail_mesh=[v for v in mesh if v>=r0**a]
            b=m.quad(parts_log,log_mesh,method='gauss-legendre',maxdegree=8)
            b+=m.quad(far_parts,tail_mesh,method='gauss-legendre',maxdegree=8)
            raw=near_direct+p0*v**(-k)/k-norm*d
            raw_parts=-(near_parts+norm*b/2)/k
            scaled=-raw/(norm*g)
            alternative=-raw_parts/(norm*g)
            discrepancy=abs(scaled-alternative)
            assert raw<0 and raw_parts<0 and discrepancy<m.mpf('1e-73')
            measured=[];remainder=scaled-1
            for j in range(1,5):
                estimate=remainder*xx**j
                measured.append(dict(order=j,exact_coefficient=str(cs[j]),prediction=m.nstr(c[j],45),
                                     extracted=m.nstr(estimate,45),difference=m.nstr(estimate-c[j],35)))
                remainder-=c[j]/xx**j
            row=dict(X=X,marchaud=m.nstr(raw,62),marchaud_parts=m.nstr(raw_parts,62),
                     normalized_negative=m.nstr(scaled,62),normalized_method_difference=m.nstr(discrepancy,35),
                     coefficients=measured,remainder_after_four=m.nstr(remainder,35),r_cutoff=m.nstr(U,35))
            rows.append(row)
        return rows


def main():
    p=argparse.ArgumentParser();p.add_argument('kappa');p.add_argument('variance');p.add_argument('out',type=Path);a=p.parse_args()
    start=time.time();cpu=time.process_time();r=evaluate(Q(a.kappa),Q(a.variance))
    a.out.write_text(json.dumps(dict(variance=a.variance,kappa=a.kappa,rows=r,pid=os.getpid(),python=platform.python_version(),
        mpmath=m.__version__,dps=100,start=start,end=time.time(),cpu_seconds=time.process_time()-cpu),indent=2)+'\n')

if __name__=='__main__':main()
