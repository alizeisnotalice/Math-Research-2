"""Exact rational certificates for a smooth negative test on a fat Cantor source."""
import argparse,json,math,os,platform,time
from fractions import Fraction as Q
from pathlib import Path


def exp_bounds(x,n):
    """Truncated exponential series and geometric tail for x>=0."""
    assert 0<=x<n+2
    terms=sum(x**k/math.factorial(k) for k in range(n+1))
    return terms,terms+(x**(n+1)/math.factorial(n+1))/(1-x/(n+2))


def atan_bounds(x,n):
    s=sum((-1)**k*x**(2*k+1)/(2*k+1) for k in range(n))
    return s,s+x**(2*n+1)/(2*n+1)


def cert():
    e=Q(1,8);s=1-e;c=1-e/2
    a=2+Q(1,5000);r=Q(3,20000);left,right=a-r,a+r
    delta=Q(1,256);h=Q(1,10**12);t=s+h
    # The right displacement is strictly inside e^(3/2)/100.
    assert left>2 and (right-2)**2<e**3/100**2 and right<Q(21,10)
    assert Q(0)<h<e and s<t<1
    mach_a=atan_bounds(Q(1,5),12);mach_b=atan_bounds(Q(1,239),12)
    pi_lo,pi_hi=16*mach_a[0]-4*mach_b[1],16*mach_a[1]-4*mach_b[0]
    assert 3<pi_lo<pi_hi<4
    assert exp_bounds(Q(1),6)[1]<3 and exp_bounds(Q(1),6)[0]>2
    # p_s/g<7/8 on I, via -log(1-e)<=e/s and exp(3e/(2s))>=1+3e/(2s).
    free_ratio=(1+3*e/(2*s))**-1
    assert free_ratio<1-e
    # Gaussian lower density and local density ratio.
    exp_cap=exp_bounds(Q(9,4),12)[1]
    assert right**2/2<Q(9,4) and exp_cap<10 and 2*pi_hi<9
    g_lower=Q(1,30)
    log_ratio=(right**2-left**2)/2
    assert log_ratio<Q(1,100) and 1/(1-log_ratio)<Q(5,4)
    # Nonnegative middle measure, using a completed-square kernel bound.
    assert e*3**2>1 and right**2/(2*s)<3 and 2*pi_lo>4
    residual_loss=Q(3)*27*Q(1,2)*(right-left)
    assert residual_loss<1-c
    # e^(-1/3)>=2/3 on the central half of the smooth bump.
    signed_margin=(c-(1-e))*Q(2,3)-c*delta*2*Q(5,4)
    assert signed_margin>Q(1,32)
    early_bound=-r*g_lower/32
    # Both component densities of Q_s are bounded above by 1/2.
    assert 2*pi_lo*s>4
    # |psi'|<4/r, and its magnitude is unimodal on either half.
    # The heat derivative is (1/2)P_h psi'', giving loss <4h/r.
    heat_loss=4*h/r
    late_bound=early_bound+heat_loss
    assert late_bound<-r/1920
    return dict(epsilon=e,s=s,c=c,center=a,radius=r,left=left,right=right,delta=delta,h=h,t=t,
       removed_length_upper=delta*(right-left),cantor_length_lower=(1-delta)*(right-left),
       residual_middle_density_lower_factor=1-c-residual_loss,free_ratio_upper=free_ratio,
       exponential_upper=exp_cap,g_density_lower=g_lower,
       normalized_early_margin=signed_margin,early_pairing_upper=early_bound,
       heat_loss_upper=heat_loss,late_pairing_upper=late_bound,simplified_late_pairing_upper=-r/1920)


def check_cantor(depth,params):
    a,b=params['left'],params['right'];delta=params['delta']
    intervals=[(a,b)]
    product=Q(1);gaps=Q(0)
    for n in range(1,depth+1):
        fraction=delta/Q(2)**n
        product*=1-fraction
        after=[]
        for lo,hi in intervals:
            half=(hi-lo)*(1-fraction)/2
            lo_gap,hi_gap=lo+half,hi-half
            assert lo<lo_gap<hi_gap<hi
            gaps+=hi_gap-lo_gap
            after.extend(((lo,lo_gap),(hi_gap,hi)))
        intervals=after
        length=sum(hi-lo for lo,hi in intervals)
        assert length==product*(b-a)==b-a-gaps
        assert len(intervals)==2**n
        assert max(hi-lo for lo,hi in intervals)<(b-a)/2**n
        assert gaps<(b-a)*delta
    return dict(stage=depth,components=len(intervals),positive_kept_length=str(length),
                exact_removed_length=str(gaps),maximum_component_length=str(max(hi-lo for lo,hi in intervals)))


def main():
    p=argparse.ArgumentParser();p.add_argument('depth',type=int);p.add_argument('out',type=Path);a=p.parse_args()
    start=time.time();cpu=time.process_time();q=cert();r=check_cantor(a.depth,q)
    rec=dict(proof={key:str(value) for key,value in q.items()},finite_cantor=r,pid=os.getpid(),
             python=platform.python_version(),start=start,end=time.time(),cpu_seconds=time.process_time()-cpu)
    a.out.write_text(json.dumps(rec,indent=2)+'\n')

if __name__=='__main__':main()
