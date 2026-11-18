"""Exact real sign check for the degree-32 shape XXBAXXBAA."""
from argparse import ArgumentParser
from fractions import Fraction
import json
from pathlib import Path
import sys

from sympy import Poly, Rational, ZZ, symbols, __version__ as sympy_version
from sympy.polys.rings import ring

sys.dont_write_bytecode=True
from census32 import descriptor, sha
from prefix_29_30 import model, slots

HERE=Path(__file__).resolve().parent
PREVIOUS=HERE.parent/'round066-gap7-next'
sys.path.insert(0,str(PREVIOUS))
from verify_c12_c16 import Interval, lower

SHAPE='XXBAXXBAA'
ROW=19


def polynomial_data(shape,row):
    parameter=descriptor(shape)
    positions=slots(32,parameter)
    R,t,*vars=ring(['t',*[f'u{k}' for k in positions]],ZZ)
    raw=dict(model(t,32,parameter,vars)[0])[row]
    T=symbols('t')
    H=Poly(T**6-4*T**5+6*T**4-7*T**3+4*T*T-2*T+1,T,domain=ZZ)
    separated={}
    for exp,c in raw.items():
        monomial=exp[1:]
        separated[monomial]=separated.get(monomial,0)+int(c)*T**exp[0]
    coefficients={m:Poly(c,T,domain=ZZ).rem(H) for m,c in separated.items()}
    coefficients={m:c for m,c in coefficients.items() if c}
    return T,H,positions,coefficients


def basis(poly):
    return [int(poly.nth(i)) for i in range(6)]


def enclosure(poly,I):
    value=Interval.exact(0)
    for c in poly.all_coeffs():
        value=value*I+int(c)
    return value


def run():
    t,H,positions,P=polynomial_data(SHAPE,ROW)
    assert positions==(8,9,12,13)
    assert set(P)=={(2,0,0,0),(1,1,0,0),(1,0,0,0),(0,2,0,0),
                    (0,1,0,0),(0,0,1,0),(0,0,0,1),(0,0,0,0)}
    bounds=(Fraction(71942965,10**8),Fraction(71942966,10**8))
    assert H.count_roots(0,1)==1 and H.count_roots(*map(Rational,bounds))==1
    I=Interval(*bounds)
    negative=[(1,1,0,0),(0,2,0,0),(0,1,0,0)]
    positive=[(0,0,1,0),(0,0,0,1)]
    assert all(enclosure(P[m],I).hi<0 for m in negative)
    assert all(enclosure(P[m],I).lo>0 for m in positive)
    a,b=lower(t)
    q=1-a[3]*b[5]
    y=1-a[4]*b[5]-t*q
    u=symbols('u')
    lower_polynomial=Poly(sum(P[m].as_expr()*u**m[0]*(y+(2*t-1)*u)**m[1]
                              for m in P if m[2]==0 and m[3]==0),u)
    assert lower_polynomial.degree()==2
    coeffs=[Poly(lower_polynomial.nth(i),t,domain=ZZ).rem(H) for i in (2,1,0)]
    boxes=[enclosure(c,I) for c in coeffs]
    for bounds,threshold in zip(boxes,(Fraction(9,10),Fraction(-9,10),Fraction(49,100)),strict=True):
        assert bounds.lo>threshold
    delta=Fraction(49,100)-Fraction(9,40)
    assert delta==Fraction(53,200)>0
    return {'status':'verified','shape':SHAPE,'row':ROW,'variable_indices':positions,
            'exact_coefficients':[[m,basis(p)] for m,p in sorted(P.items(),reverse=True)],
            'post_substitution_basis':[basis(c) for c in coeffs],
            'post_substitution_boxes':[[str(b.lo),str(b.hi)] for b in boxes],
            'strict_lower_bound':str(delta),'python':sys.version,'sympy':sympy_version,
            'sources_sha256':{'verify_real32.py':sha(Path(__file__)),
                              'verify_c12_c16.py':sha(PREVIOUS/'verify_c12_c16.py')}}


def main():
    if not __debug__:
        raise SystemExit('assertions required')
    p=ArgumentParser(description=__doc__)
    p.add_argument('--obstruction',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    r=run()
    q=json.loads(a.obstruction.read_bytes())
    assert q['shape']==r['shape'] and list(r['variable_indices'])==q['exact_system']['variable_indices']
    original=dict(q['exact_system']['early_zero_equations'])[r['row']]
    assert json.loads(json.dumps(r['exact_coefficients']))==original
    r['obstruction_sha256']=sha(a.obstruction)
    with a.out.open('x') as f:
        f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(r['shape'],r['row'],r['strict_lower_bound'],flush=True)


if __name__=='__main__':
    main()
