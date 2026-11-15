#!/usr/bin/env python3
"""Independent exact scalar certificates for the two R38 Green-domain proofs."""
from fractions import Fraction as F
from pathlib import Path
import datetime
import hashlib
import json
import os
import platform
import subprocess
import sys


def fs(x):
    return str(x.numerator)+'/'+str(x.denominator)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def add(*polys):
    out={}
    for poly in polys:
        for k,v in poly.items():
            out[k]=out.get(k,F(0))+v
    return {k:v for k,v in out.items() if v}


def mul(a,b):
    out={}
    for i,x in a.items():
        for j,y in b.items():
            out[i+j]=out.get(i+j,F(0))+x*y
    return {k:v for k,v in out.items() if v}


def scale(poly,c):
    return {k:v*c for k,v in poly.items() if v*c}


def deriv(poly):
    return {k-1:k*v for k,v in poly.items() if k}


def atan_bounds(x,n):
    assert F(0)<x<F(1) and n>=1
    total=sum((F(-1) if k%2 else F(1))*x**(2*k+1)/F(2*k+1)
              for k in range(n))
    next_abs=x**(2*n+1)/F(2*n+1)
    return (total,total+next_abs) if n%2==0 else (total-next_abs,total)


def qadd(x,y):
    return (x[0]+y[0],x[1]+y[1])


def qneg(x):
    return (-x[0],-x[1])


def qmul(x,y):
    return (x[0]*y[0]+2*x[1]*y[1],x[0]*y[1]+x[1]*y[0])


def qdiv(x,y):
    den=y[0]*y[0]-2*y[1]*y[1]
    assert den
    return qmul(x,(y[0]/den,-y[1]/den))


def main():
    assert os.path.ismount('/Volumes/data'), 'data volume not mounted'
    base=Path(__file__).resolve().parent.parent
    author=Path('/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-astra-kernel/planar-centered-square-weak11/proofs')
    source_hashes={
        'killed-green-dihedral-domain-obstruction.md':'16718b81bdcd80af4ea0b6d2bfb261ee2e586de08f45e19eedaf5580d015314a',
        'killed-green-homothetic-diamond-obstruction.md':'1bc3087156ef3f34f03a0e5734fb7a9b3e04aef33f0603821dccae389ca25d05'}
    # Hash checks bind this review to the exact submitted author versions.
    for name,digest in source_hashes.items():
        assert sha(author/name)==digest,(name,sha(author/name))

    lhs={4:F(1),2:F(1),1:F(-12,5),0:F(1)}
    rhs=add(mul({2:F(1),0:F(-1,3)},{2:F(1),0:F(-1,3)}),
            scale(mul({1:F(1),0:F(-18,25)},{1:F(1),0:F(-18,25)}),F(5,3)),
            {0:F(28,1125)})
    assert lhs==rhs
    assert F(28,1125)>0

    # Quotient derivative for h(t)=t(1-t)/(1+t).
    num={1:F(1),2:F(-1)}; den={0:F(1),1:F(1)}
    h_derivative_num=add(mul(deriv(num),den),scale(mul(num,deriv(den)),F(-1)))
    assert h_derivative_num=={0:F(1),1:F(-2),2:F(-1)}
    one=(F(1),F(0)); tstar=(F(-1),F(1))
    critical=qadd(qadd(qmul(tstar,tstar),qmul((F(2),F(0)),tstar)),(F(-1),F(0)))
    assert critical==(F(0),F(0))
    h_at_star=qdiv(qmul(tstar,qadd(one,qneg(tstar))),qadd(one,tstar))
    assert h_at_star==qmul(tstar,tstar)==(F(3),F(-2))

    # Machin identity: (5+i)^4(239-i)=114244(1+i).
    # Its angle lies in (0,pi/2), hence pi=16atan(1/5)-4atan(1/239).
    z=(1,0)
    for _ in range(4):
        z=(z[0]*5-z[1],z[0]+z[1]*5)
    z=(239*z[0]+z[1],239*z[1]-z[0])
    assert z==(114244,114244)
    lo5,hi5=atan_bounds(F(1,5),10)
    lo239,hi239=atan_bounds(F(1,239),4)
    pi_lo=16*lo5-4*hi239
    pi_hi=16*hi5-4*lo239
    assert 0<pi_lo<pi_hi
    assert pi_lo>F(31415,10000)>3
    sqrt_lo=F(14142,10000); sqrt_hi=F(14143,10000)
    assert sqrt_lo*sqrt_lo<2<sqrt_hi*sqrt_hi
    d4_lo=F(6,5)*pi_lo
    hom_lo=pi_lo*(1+sqrt_lo)/2
    simple_hom_lo=F(31415,10000)*(1+F(14142,10000))/2
    assert simple_hom_lo==F(75842093,20000000)
    assert hom_lo>simple_hom_lo>F(379,100)
    assert d4_lo>F(18,5)>F(3590865312933,10**12)

    parseval_rows=[]
    for n in range(8):
        k=4*n+1
        radial_log_integral=F(1,(2*k)**2)
        coefficient=k*k*radial_log_integral
        assert coefficient==F(1,4)
        parseval_rows.append({'n':n,'k':k,'radial_log_integral':fs(radial_log_integral),
                              'k_squared_times_integral':fs(coefficient)})
    gitbase=subprocess.run(['git','-C',str(base.parent),'rev-parse','HEAD'],
                           text=True,capture_output=True,check=True).stdout.strip()
    review=base/'reviews/killed-green-domain-second-review.md'
    out={'status':'PASS: exact scalar certificates; analytic review supplies geometric claims',
         'timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
         'python':sys.version,'executable':sys.executable,'platform':platform.platform(),
         'arithmetic':'fractions.Fraction and Q(sqrt(2)), exact','workers':1,
         'random_seed':None,'floating_tolerance':None,'git_base':gitbase,
         'script_sha256':sha(Path(__file__).resolve()),
         'author_source_sha256':source_hashes,
         'review_sha256':sha(review) if review.exists() else None,
         'sos_coefficients':{str(k):fs(v) for k,v in sorted(rhs.items())},
         'sos_positive_remainder':'28/1125',
         'h_derivative_numerator':{str(k):fs(v) for k,v in sorted(h_derivative_num.items())},
         'h_max_exact_Qsqrt2':{'rational_part':'3/1','sqrt2_part':'-2/1'},
         'machin_complex_product':z,'atan_terms':{'one_fifth':10,'one_239th':4},
         'pi_rational_bounds':{'lower':fs(pi_lo),'upper':fs(pi_hi)},
         'sqrt2_rational_bounds':{'lower':fs(sqrt_lo),'upper':fs(sqrt_hi)},
         'd4_constant_lower':fs(d4_lo),'homothetic_constant_lower':fs(hom_lo),
         'quoted_simple_homothetic_lower':fs(simple_hom_lo),
         'parseval_normalization':parseval_rows,
         'scope':'No domain sampling, no conformal numerical approximation, no non-D4 group averaging, no Green minimum computation, no R33 or R9 rerun. Exact arithmetic checks do not replace Riemann-map, a.e.-vertex, or mixture proofs.'}
    path=base/'computation/results/killed-green-r038-second-review/exact-scalars.json'
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':out['status'],'output':str(path),'pi_gt_31415_over_10000':True,
                      'D4_constant_gt_18_over_5':True,'homothetic_constant_gt_379_over_100':True,
                      'source_hashes_verified':True},indent=2))

if __name__=='__main__':
    main()
