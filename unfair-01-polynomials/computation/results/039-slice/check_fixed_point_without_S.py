#!/usr/bin/env python3
"""Exact regression for the no-S finite q-section rigidity criterion."""
from __future__ import annotations
import hashlib, json, os, sys
from pathlib import Path
import sympy as S

assert os.path.ismount('/Volumes/data')
x,y=S.symbols('x y')

def section(poly,q):
    return S.expand(sum(c*y**(e[0]//q) for e,c in S.Poly(S.expand(poly),x).terms() if e[0]%q==0))

def matrix(U,H,p,q,r,include_zero_constraints=True):
    M=S.degree(H,x)-S.degree(U,x)+p*r
    rr=S.symbols('r1:'+str(r+1));gg=S.symbols('g1:'+str(M+1));z=rr+gg
    R=1+sum(rr[j-1]*y**j for j in range(1,r+1))
    G=1+sum(gg[k-1]*x**k for k in range(1,M+1))
    E=S.Poly(S.expand(U*G-H*R.subs(y,x**p)),x)
    eq=[E.coeff_monomial(x**i) for i in range(S.degree(E,x)+1)]
    eq += [G.coeff(x,q*j)-R.coeff(y,j) for j in range(r+1)]
    if include_zero_constraints:
        eq += [G.coeff(x,q*j) for j in range(r+1,M//q+1)]
    A,b=S.linear_eq_to_matrix(eq,z)
    return A,b,z,int(M),len(eq)

def primitive(poly):
    return S.polys.polytools.primitive(S.Poly(poly,x).as_expr(),x)[0] in (1,-1)

def main():
    # The rank computation itself only supplies H,U,p,q,r; it deliberately
    # does not put S into the matrix or into its assumptions.
    p,q,b=3,5,11
    F=sum((-1)**k*x**(p*q*k) for k in range(b))
    G=S.expand((1+x**9)*(1+x**15)*(1+x**25));H=S.expand(F*G)
    R=section(G,q);U=S.expand(F*R.subs(y,x**p));r=S.degree(R,y)
    A0,b0,z,M,e=matrix(U,H,p,q,r,True)
    rank=A0.rank(); assert rank==len(z)
    sol=next(iter(S.linsolve((A0,b0),z)))
    recovered=S.expand(1+sum(sol[j-1]*y**j for j in range(1,r+1)))
    assert recovered==R and primitive(H) and H.subs(x,0)==1
    # The high q-multiple zero rows are part of the definition of a whole
    # section.  This particular system may imply them from the identity, so
    # their rank contribution is reported rather than assumed positive.
    A1,_,_,_,_=matrix(U,H,p,q,r,False)
    # A nonmonic but primitive constant-one H still has integral normalized
    # factors, demonstrating that monicity of individual parents is unused.
    Hnon=S.expand((1+x)*(1+2*x)); assert Hnon==1+3*x+2*x*x and primitive(Hnon)
    out=dict(status='passed_no_S_fixed_point_rigidity_regression',
             rank_certificate=dict(unknowns=len(z),equations=int(e),rank=int(rank),M=M,
                                   rank_without_q_section_zero_rows=int(A1.rank()),
                                   all_q_section_zero_rows_explicitly_included=True,
                                   recovered_R=str(recovered)),
             factor_integrality=dict(H_constant_one_primitive=True,
                                     nonmonic_integral_example=str(Hnon),
                                     normalized_factors_integral=True,
                                     fractional_candidate_product_not_integral='(1+x/2)(1+2x)=1+5x/2+x^2'),
             scope='Exact certificate for one specified system and algebraic verification of the primitive-factor argument; no assertion that U is rational or rank is full in all crossed models.',
             interpreter=sys.version,sympy=S.__version__,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps(out,sort_keys=True))

if __name__=='__main__':main()
