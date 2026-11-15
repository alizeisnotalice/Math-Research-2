#!/usr/bin/env python3
"""Exact symbolic interfaces and independent 1D quadratures for Cauchy atom kernels."""
from pathlib import Path
from fractions import Fraction
import datetime
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
import sympy as sp
import mpmath as mp


def main():
    assert os.path.ismount('/Volumes/data'), 'data volume not mounted'
    started=time.monotonic()
    base=Path(__file__).resolve().parent.parent
    r,d,R,k=sp.symbols('r d R k',positive=True)
    q=sp.symbols('theta',real=True)
    checks={}
    def verify(name,expr):
        result=sp.factor(sp.together(expr))
        assert result==0,(name,result)
        checks[name]='exact symbolic zero'
    w4=1/r-R/(2*(R**2-d**2))-1/(2*(R+r))
    ray4=-1/r**2-R*k**2*r/(R**2-k**2*r**2)**2+1/(2*(R+r)**2)
    verify('four_ray_derivative',sp.diff(w4.subs(d,k*r),r)-ray4)
    m4=1-(2*R-1)/(2*(R**2-1))
    verify('four_vertex_value',w4.subs({r:1,d:1})-m4)
    f4=R/m4
    g4=2*R**4-4*R**3-R**2+1
    verify('four_radius_derivative',sp.diff(f4,R)-2*g4/(2*R**2-2*R-1)**2)
    wa=1/r-1/(2*(1-d**2))-1/(2*(1+r))
    verify('core_sign_numerator',wa-(2-r**2-(r+2)*d**2)/(2*r*(1+r)*(1-d**2)))
    root0=sp.sqrt(3)-1
    verify('core_crossover',root0**3+3*root0**2-2)
    wb=1/r-1/(4*(r-1))-1/(2*(r+1))-1/(4*(d+1))
    den=r**2+r-4
    p=den/(4*r*(r**2-1))
    threshold=r-2+4*(r-1)/den
    verify('tail_partial_fraction',wb-(p-1/(4*(d+1))))
    verify('tail_zero_boundary',wb.subs(d,threshold))
    z=2*(r-2)/(r*(r-1))
    verify('tail_width',r-threshold-2*(r-2)*(r+1)/den)
    verify('tail_integral_linear_piece',4*p*(r-threshold)-z)
    verify('tail_integral_log_argument',(r+1)/(threshold+1)-(1+z))
    a=sp.Rational(5,2)-sp.sqrt(17)/2
    b=sp.Rational(5,2)+sp.sqrt(17)/2
    verify('tail_root_product',a*b-2)
    verify('tail_root_sum',a+b-5)

    wc=1/r-R/(2*(R**2-d**2))-q/(2*(R+r))-(1-q)*R/(2*(R**2-r**2))
    wo=(3+q)/(4*r)-q*R/(2*(R**2-d**2))-(1-q)/(4*d)-(1-q)/(4*(R+d))-(1+q)/(4*(R+r))
    verify('eight_branch_continuity',wc.subs(d,R-r)-wo.subs(d,R-r))
    second=-q*R*(R**2+3*d**2)/(R**2-d**2)**3-(1-q)/(2*d**3)-(1-q)/(2*(R+d)**3)
    verify('eight_outer_d_concavity',sp.diff(wo,d,2)-second)
    cross=3/(4*r)-1/(4*(2*R-r))-q/(2*(R+r))-(1-q)*R/(2*(R**2-r**2))
    verify('eight_cross_line_value',wc.subs(d,R-r)-cross)
    cross_derivative=-3/(4*r**2)-1/(4*(2*R-r)**2)+q/(2*(R+r)**2)-(1-q)*R*r/(R**2-r**2)**2
    verify('eight_cross_line_derivative',sp.diff(cross,r)-cross_derivative)
    axis_out=(1+q)/(2*r)-q*R/(2*(R**2-r**2))-1/(2*(R+r))
    verify('eight_outer_axis_value',wo.subs(d,r)-axis_out)
    axis_out_derivative=-(1+q)/(2*r**2)-q*R*r/(R**2-r**2)**2+1/(2*(R+r)**2)
    verify('eight_outer_axis_derivative',sp.diff(axis_out,r)-axis_out_derivative)
    axis_in=1/r-(2-q)*R/(2*(R**2-r**2))-q/(2*(R+r))
    verify('eight_inner_axis_value',wc.subs(d,r)-axis_in)
    axis_in_derivative=-1/r**2-(2-q)*R*r/(R**2-r**2)**2+q/(2*(R+r)**2)
    verify('eight_inner_axis_derivative',sp.diff(axis_in,r)-axis_in_derivative)
    A=q*m4+(1-q)*R/(2*(R+1))
    C=1-R/(2*(2*R-1))-q/(2*(R+1))-(1-q)*R/(2*(R**2-1))
    verify('eight_axis_candidate',axis_out.subs(r,1)-A)
    verify('eight_cross_candidate',cross.subs(r,1)-C)
    mlarge=1-(2*R-q)/(2*(R**2-1))
    verify('eight_large_radius_minimum',wc.subs({r:1,d:1})-mlarge)
    theta_c=(R**2+R-1)/(2*R**2+R-1)
    difference=(2-R)*(R**2+R-1-q*(2*R**2+R-1))/(2*(R**2-1)*(2*R-1))
    verify('eight_switch_weight',A-C-difference)
    g8=2*R**4-4*R**3+(3*q-4)*R**2+2-q
    verify('eight_large_radius_derivative',sp.diff(R/mlarge,R)-2*g8/(2*R**2-2*R+q-2)**2)
    switch_gap=sp.factor(A.subs(q,theta_c)/R-(2+theta_c)/12)
    # Purely analytic four-axis exclusion, proposed by Astra and audited here.
    x=sp.symbols('x',real=True)
    P=10*R**3-38*R**2+28*R+19
    verify('four_cost_factor_polynomial',R/(2*m4)-sp.Rational(19,10)-P/(20*m4*(R**2-1)))
    shifted=10*x**3+22*x**2-4*x+3
    verify('four_cost_factor_shift',P.subs(R,x+2)-shifted)
    verify('four_cost_factor_positive_halfline',shifted-(10*x**3+22*(x-sp.Rational(1,11))**2+sp.Rational(31,11)))
    first_piece=1+sp.log(3)/2
    second_piece=2*sp.log(2)-sp.log(3)/2-sp.expand_log(sp.log(sp.Rational(4,3)),force=True)
    verify('truncated_core_integral',first_piece+second_piece-1-sp.log(3))
    log3_lower=Fraction(1)+Fraction(1,12)+Fraction(1,80)
    exact_cost_floor=Fraction(19,10)*(1+log3_lower)
    verify('eight_small_radius_axis_cap_four', (R-1)/R-m4-(2-R)/(2*R*(R**2-1)))
    verify('eight_small_radius_axis_cap_amgm', R/4-(R-1)/R-(R-2)**2/(4*R))
    verify('eight_small_radius_diag_cap', R/4-R/(2*(R+1))-R*(R-1)/(4*(R+1)))
    t=sp.symbols('t',positive=True)
    verify('rectangle_potential_primitive',sp.diff(t*sp.log(t)-t,t,2)-1/t)
    assert log3_lower==Fraction(263,240)
    assert exact_cost_floor==Fraction(9557,2400)>Fraction(199,50)

    def stable_z_log(z):
        if abs(z)>mp.mpf('0.01'):
            return z-mp.log1p(z)
        total=mp.mpf(0); power=z*z; n=2
        while True:
            term=power/n if n%2==0 else -power/n
            total+=term
            if abs(term)<=mp.eps*max(abs(total),mp.mpf('1e-1000')):
                return total
            power*=z;n+=1
            assert n<10000

    numeric=[]
    for precision in [70,110]:
        mp.mp.dps=precision
        sq2=mp.sqrt(2); r0=mp.sqrt(3)-1
        def core_g(rr,D):
            if rr==0:return mp.mpf(4)
            return 4*(D/rr-mp.atanh(D)/2-D/(2*(1+rr)))
        c1=mp.quad(lambda rr:core_g(rr,rr),[0,r0])
        tmax=mp.sqrt(sq2-r0)
        def transformed(t):
            rr=sq2-t*t
            cap=mp.sqrt(max(mp.mpf(0),(2-rr*rr)/(rr+2)))
            return 2*t*core_g(rr,cap)
        c2=mp.quad(transformed,[0,tmax])
        core_r=c1+c2
        def core_d_transform(t):
            if t==0:return mp.mpf(0)
            dd=r0*t*t
            rho=(-dd*dd+mp.sqrt(dd**4-8*dd*dd+8))/2
            inner=mp.log(rho/dd)-mp.log((1+rho)/(1+dd))/2-(rho-dd)/(2*(1-dd*dd))
            return 8*r0*t*inner
        core_d=mp.quad(core_d_transform,[0,1])
        aa=(5-mp.sqrt(17))/2;bb=(5+mp.sqrt(17))/2
        tail_closed=aa*mp.log(aa)+bb*mp.log(bb)-6*mp.log(2)-2
        def tail_inverse(t):
            if t==0:return mp.mpf(2)
            zz=2*t*(1-2*t)/(1-t)
            return stable_z_log(zz)/(t*t)
        tail_quad=mp.quad(tail_inverse,[0,mp.mpf(1)/4,mp.mpf(1)/2])
        I=core_r+tail_closed
        optR=mp.findroot(lambda rr:2*rr**4-4*rr**3-rr**2+1,(2,mp.mpf('2.25')))
        m=1-(2*optR-1)/(2*(optR**2-1))
        J=optR*I/(2*m)
        assert abs(core_r-core_d)<mp.mpf(10)**(-(precision-15))
        assert abs(tail_quad-tail_closed)<mp.mpf(10)**(-(precision-15))
        numeric.append({'digits':precision,'core_by_r':mp.nstr(core_r,precision),
                        'core_by_d':mp.nstr(core_d,precision),
                        'tail_closed':mp.nstr(tail_closed,precision),
                        'tail_quadrature':mp.nstr(tail_quad,precision),
                        'I_four':mp.nstr(I,precision),'optimal_R':mp.nstr(optR,precision),
                        'm_opt':mp.nstr(m,precision),'cost_J_four_numeric':mp.nstr(J,precision),
                        'core_discrepancy':mp.nstr(abs(core_r-core_d),8),
                        'tail_discrepancy':mp.nstr(abs(tail_quad-tail_closed),8)})
    # Exact selected rational checks against a direct eight-charge evaluator.
    def direct_w(x,y,RR,theta):
        def u(a,b):return Fraction(1,1)/(abs(a)+abs(b))
        axis=[(RR,0),(-RR,0),(0,RR),(0,-RR)]
        B=RR/2
        diag=[(B,B),(B,-B),(-B,B),(-B,-B)]
        return u(x,y)-theta*sum(u(x-a,y-b) for a,b in axis)/4-(1-theta)*sum(u(x-a,y-b) for a,b in diag)/4
    rational=[]
    for RR in [Fraction(5,4),Fraction(3,2),Fraction(2),Fraction(5,2)]:
        for theta in [Fraction(0),Fraction(1,2),Fraction(1)]:
            av=theta*(1-(2*RR-1)/(2*(RR**2-1)))+(1-theta)*RR/(2*(RR+1)) if RR<2 else 1-(2*RR-theta)/(2*(RR**2-1))
            got=direct_w(Fraction(1),Fraction(0),RR,theta)
            assert got==av
            if RR<2:
                cv=1-RR/(2*(2*RR-1))-theta/(2*(RR+1))-(1-theta)*RR/(2*(RR**2-1))
                gotcross=direct_w(RR/2,1-RR/2,RR,theta)
                assert gotcross==cv
            else:gotcross=None
            rational.append({'R':str(RR),'theta':str(theta),'axis_value':str(got),
                             'cross_value':None if gotcross is None else str(gotcross)})
    gitbase=subprocess.run(['git','-C',str(base.parent),'rev-parse','HEAD'],capture_output=True,text=True,check=True).stdout.strip()
    out={'status':'PASS exact algebra and actual 1D numerical cross-checks; numerical values are not directed certificates',
         'timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
         'python':sys.version,'executable':sys.executable,'platform':platform.platform(),
         'sympy':sp.__version__,'mpmath':mp.__version__,'workers':1,'git_base':gitbase,
         'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'symbolic_checks':checks,'eight_switch_gap_over_radius':str(switch_gap),
         'rational_checks':rational,'numerical_cross_checks':numeric,
         'exact_four_axis_exclusion':{'I_strictly_above':'1+log(3)',
              'R_over_two_m_strictly_above':'19/10','log3_strictly_above':str(log3_lower),
              'J_strictly_above':str(exact_cost_floor),'strictly_above_3_98':True},
         'elapsed_seconds':time.monotonic()-started,
         'scope':'Analytic proof establishes geometry and exact minimum formulas. CAS checks are scalar identities only. mpmath is nondirected numerical evidence. No large 2D integration, MATLAB, LP, or frozen proof rerun.'}
    path=base/'computation/results/cauchy-positive-part-atom-interfaces/checks.json'
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':out['status'],'symbolic_check_count':len(checks),
                      'eight_switch_gap_over_radius':str(switch_gap),
                      'last_numerical':numeric[-1],'elapsed_seconds':out['elapsed_seconds'],
                      'output':str(path)},indent=2))

if __name__=='__main__':main()
