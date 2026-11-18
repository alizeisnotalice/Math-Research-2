"""Exactly exclude the degree-31 shape XXBAABBB by its 23rd-row square."""
from argparse import ArgumentParser
from fractions import Fraction
from hashlib import file_digest
import json
from pathlib import Path
import sys

from sympy import Poly, Rational, ZZ, symbols, __version__ as sympy_version
from verify_c12_c16 import Interval, lower

SHAPE='XXBAABBB'
RULES='XXBAABBBBBBAABZZ'


def vector(poly):
    return [int(poly.nth(k)) for k in range(6)]


def evaluate(coefficients,t):
    answer=0
    for x in coefficients:
        answer=answer*t+int(x)
    return answer


def run():
    t,u,v=symbols('t u v')
    H=Poly(t**6-4*t**5+6*t**4-7*t**3+4*t*t-2*t+1,t,domain=ZZ)
    assert RULES==SHAPE+''.join({'A':'A','B':'B','X':'Z'}[x] for x in SHAPE[::-1])
    A,B=lower(t)
    unknown=iter((u,v))
    for k,state in enumerate(RULES,8):
        inside=sum(a*B[k-i] for i,a in enumerate(A) if 1<=k-i<len(B) and i>0)
        S=1-inside
        if state=='A':
            x,y=S,0
        elif state=='B':
            x,y=0,S
        elif state=='X':
            x=next(unknown)
            y=S-x
        else:
            assert state=='Z'
            x=y=0
        A.append(x)
        B.append(y)
    assert k==23 and state=='Z' and not list(unknown)
    F=Poly(S,u,v)
    expected_monomials=(u**2,u*v,u,v**2,v,1)
    numbers=[Poly(F.coeff_monomial(mon),t,domain=ZZ).rem(H) for mon in expected_monomials]
    assert F.total_degree()==2 and sum(bool(g) for g in numbers)==len(F.terms())
    a,b,d,c,e,f=numbers
    assert f==Poly(1,t,domain=ZZ)
    assert (4*c*a-b*b).rem(H).is_zero
    assert (4*c*d-2*b*e+4*c).rem(H).is_zero
    t_bounds=(Fraction(71942965,10**8),Fraction(71942966,10**8))
    assert H.count_roots(0,1)==1 and H.count_roots(*map(Rational,t_bounds))==1
    t_interval=Interval(*t_bounds)
    lowA,lowB=lower(t_interval)
    q=1-lowA[3]*lowB[5]
    c_box=evaluate(c.all_coeffs(),t_interval)
    e_box=evaluate(e.all_coeffs(),t_interval)
    assert 0<q.lo<=q.hi<Fraction(9,10)
    assert c_box.lo>2
    assert 0<e_box.lo<=e_box.hi<Fraction(1,2)
    bound=1-Fraction(9,10)-Fraction(1,32)
    assert bound==Fraction(11,160)>0
    return {'status':'verified','shape':SHAPE,'row23_coefficients':[vector(x) for x in numbers],
            'monomials':[str(mon) for mon in expected_monomials],
            'rational_boxes':{k:[str(x.lo),str(x.hi)] for k,x in (('q',q),('c',c_box),('e',e_box))},
            'strict_lower_bound':str(bound),'python':sys.version,'sympy':sympy_version}


def main():
    if not __debug__:
        raise SystemExit('assertions required')
    p=ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True)
    args=p.parse_args()
    data=run()
    sources={}
    here=Path(__file__).resolve().parent
    for n in ('verify_real31.py','verify_c12_c16.py'):
        with (here/n).open('rb') as f:
            sources[n]=file_digest(f,'sha256').hexdigest()
    data['sources_sha256']=sources
    with args.out.open('x') as out:
        out.write(json.dumps(data,sort_keys=True,indent=2)+'\n')
    print(data['shape'],data['strict_lower_bound'],flush=True)


if __name__=='__main__':
    main()
