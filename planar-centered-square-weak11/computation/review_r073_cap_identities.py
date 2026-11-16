#!/usr/bin/env python3
"""Independent exact rational-identity checks for the new cap construction."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys
import time
import sympy as s


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();start=time.perf_counter()
    n,k,l,d,R=s.symbols('n k l d R',positive=True);q=s.symbols('q',nonnegative=True)
    t=n/(n+1);c=t/(1+2*d);ak=((n+1)*k-n)/n**2;al=((n+1)*l-n)/n**2;P=(1-ak)*(1-al);G=(n-k)*(n-l)/n**3
    E=(R+2*d)**2-(ak+2*d)*(al+2*d);f=k*l/n**2-t*ak*al
    identities={
      'expanded_overlap_compensation':E-((1+2*d)*(R**2-ak*al)+2*d*(P-(R-1)**2)),
      'cap_gap_equals_nP':1-t-f-n*P/(n+1)**2,
      'full_excess_factorization':k*l/n**2+c*E-R**2-(-(1-2*n*d)*G/(1+2*d)-(R-1)*((1-c)*(R+1)-4*c*d)),
      'second_factor_positive_form':(1-c)*(R+1)-4*c*d-((1-c)*(R-1)+2*(1-t)),
      'zero_capture_density':c*(1+2*d)**2-t*(1+2*d),
      'outer_singleton_mass_above_5_over_2':t*(3-6/n**2)-s.Rational(5,2)-(n**2-5*n-12)/(2*n*(n+1))}
    checks={}
    for name,value in identities.items():
        numerator,denominator=s.fraction(s.cancel(value));assert s.Poly(numerator,n,k,l,d,R).is_zero;checks[name]=dict(exact_numerator=str(numerator),denominator=str(denominator),passed=True)
    pair=s.Poly(s.expand((n*n-4*n-4).subs(n,q+5)),q);outer=s.Poly(s.expand((n*n-5*n-12).subs(n,q+7)),q)
    assert pair.all_coeffs()==[1,6,1] and outer.all_coeffs()==[1,9,2]
    nn=9;mesh=100;eps=F(1,nn*nn);delta=F(1,mesh);density=F(nn,nn+1)/(1+2*delta);lam=F(9999,10000);atom=density*delta**2*lam
    assert delta<=F(1,2*nn) and density==F(15,17) and atom==F(29997,340000000)
    outside=300**2-105**2;local=198**2-105**2;full_lower=(1+eps)**2+local*atom;target=(2-eps)**2-F(7,2);outer_margin=local*atom-target
    assert outside==78975 and outside+nn*nn==79056 and local==28179
    assert full_lower==F(7832077922743,2230740000000)>F(7,2) and outer_margin>0
    total=1+outside*atom;gap=F(299,100)**2-F(35,33)*total;assert total==F(108360523,13600000) and gap==F(8787943,17952000)>0
    result=dict(status='PASS_INDEPENDENT_EXACT_CAP_RATIONAL_IDENTITIES',source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),python=sys.version,sympy=s.__version__,reviewed_domains='n integer>=2;1<=k,l<=n;0<delta<=1/(2n);R>=1 for cap factorization;zero-capture andR<=1 separately proved',rational_identities=checks,strict_band_gap_n_ge5_positive_polynomial=str(pair.as_expr()),outer_residual_n_ge7_positive_polynomial=str(outer.as_expr()),explicit_n9=dict(outside_atoms=outside,total_atoms=outside+81,local_complete_cells=local,mass_per_outside_atom=str(atom),local_full_halo_area_lower=str(full_lower),local_full_halo_above_7_over_2=str(full_lower-F(7,2)),strong_outer_gap_already_paid_by_residual_singletons=True,outer_residual_singleton_credit_minus_target=str(outer_margin),nonterminal_gap=str(gap)),elapsed_seconds=time.perf_counter()-start,scope='Exact formal rational identity numerators and positive polynomial coefficients; signs/domain and source geometry separately proved in review. Explicit finite input given by counts/formula, no79056source/subset/halo enumeration. Localfull obstruction does not violate actual deletion or outer-gap target.')
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(status=result['status'],identities=len(checks),n9_lower=str(full_lower),seconds=result['elapsed_seconds'])))


if __name__=='__main__':main()
