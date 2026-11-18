"""High-precision quadrature for the reverse-wing expansion on 0 < x <= e^-2."""
import argparse,json,os,platform,time
from pathlib import Path
import mpmath as m


def evaluate(Li):
    with m.workdps(95):
        L=m.mpf(Li); b=m.mpf('0.3'); d=m.mpf('0.2'); e2=m.exp(-2*L)
        def q(u):
            ell=L-m.log(u)
            return u*u/4/(ell+b*m.log(ell)+d)+u*u*e2/2
        def direct(u):
            if not u:
                return m.mpf(1)
            return m.exp(-u-q(u))
        def difference(u):
            return m.exp(-u)*(m.expm1(-q(u))+u*u/(4*L))
        knots=[m.mpf(x) for x in ('0','0.125','0.5','1','2','4','8','16','32','64','128','224')]
        val=m.quad(direct,knots,method='tanh-sinh',maxdegree=9)
        r_direct=L*L*(val-1+1/(2*L))
        # The logarithmic coordinate removes the endpoint log from the independent quadrature.
        left=[m.log(8),4,8,16,32,64,128,240]
        r_diff=L*L*(m.quad(lambda v:m.exp(-v)*difference(m.exp(-v)),left,
                            method='gauss-legendre',maxdegree=8)
                   +m.quad(difference,knots[1:]+[256],method='gauss-legendre',maxdegree=8))
        target=(b*m.log(L)+d+m.euler)/2
        A=b*m.log(L)+d
        psi3,psi5=m.digamma(3),m.digamma(5)
        third=-((A-psi3)**2+m.polygamma(1,3)+b*psi3)/2-3*(A-psi5)/2-m.mpf(15)/8
        abs_difference=abs(r_direct-r_diff)
        assert abs_difference<m.mpf('1e-76')
        tail_direct=L*L*m.exp(-224)
        tail_diff=L*L*m.exp(-256)*(1+(256**2+2*256+2)/(4*L))
        tail_diff+=L*L*(1/(2*L)+e2/2)*m.exp(-720)/3
        fields=dict(L=Li,b=str(b),d=str(d),work_dps=m.mp.dps,direct_second_remainder=r_direct,
                    cancellation_aware_second_remainder=r_diff,predicted_second_remainder=target,
                    second_remainder_minus_prediction=r_diff-target,relative_difference=(r_diff-target)/target,
                    quadrature_absolute_difference=abs_difference,direct_truncation_bound=tail_direct,
                    stable_truncation_bound=tail_diff,predicted_third_coefficient=third,
                    second_remainder_corrected_by_third=r_diff-target-third/L)
        return {k:m.nstr(v,64) if isinstance(v,m.mpf) else v for k,v in fields.items()}


def main():
    p=argparse.ArgumentParser();p.add_argument('L',type=int);p.add_argument('out',type=Path);a=p.parse_args()
    start=time.time();cpu=time.process_time();result=evaluate(a.L)
    result.update(pid=os.getpid(),start_epoch=start,end_epoch=time.time(),cpu_seconds=time.process_time()-cpu,
                  python=platform.python_version(),mpmath=m.__version__)
    a.out.write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
