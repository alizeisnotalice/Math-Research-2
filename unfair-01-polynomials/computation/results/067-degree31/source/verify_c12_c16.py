"""Exact independent check that the seventh gap forces C12 through C16 equal one."""
from argparse import ArgumentParser
from dataclasses import dataclass
from fractions import Fraction
from hashlib import file_digest
import json
from pathlib import Path
import sys

from sympy import Poly, Rational, ZZ, __version__ as sympy_version, symbols


@dataclass(frozen=True)
class Interval:
    lo: Fraction
    hi: Fraction

    @staticmethod
    def exact(value):
        if isinstance(value, Interval):
            return value
        return Interval(Fraction(value), Fraction(value))

    def __add__(self, right):
        other=self.exact(right)
        return Interval(self.lo+other.lo,self.hi+other.hi)

    __radd__=__add__

    def __neg__(self):
        return Interval(-self.hi,-self.lo)

    def __sub__(self,right):
        return self+-self.exact(right)

    def __rsub__(self,left):
        return self.exact(left)+-self

    def __mul__(self,right):
        other=self.exact(right)
        products=[x*y for x in (self.lo,self.hi) for y in (other.lo,other.hi)]
        return Interval(min(products),max(products))

    __rmul__=__mul__


A_KNOWN=(0,1,3,4)
B_KNOWN=(0,1,2,5)
FORCED={12:((10,11,12),(8,9,11,12)),
        13:((8,11,12,13),(9,10,12,13)),
        14:((9,12,13,14),(10,11,13,14)),
        15:((10,13,14,15),(11,12,14,15)),
        16:((11,14,15,16),(12,13,15,16))}


def lower(t):
    m=1-t+t*t
    n=1-t*m
    p=1-(1-t)*n
    v=t*t*m
    return [1,t,0,n,p,0,0,0],[1,1-t,m,0,0,v,0,0]


def implied_zeros(k):
    a=tuple(sorted(k-i for i in B_KNOWN if k-i>=8))
    b=tuple(sorted(k-i for i in A_KNOWN if k-i>=8))
    return a,b


def forced_previous(k,t,free):
    A,B=lower(t)
    deadA,deadB=implied_zeros(k)
    choices=iter(free)
    for j in range(8,k):
        internal=sum(A[i]*B[j-i] for i in range(1,j))
        if j in deadA and j in deadB:
            assert j==k-1 and not list(choices)
            return internal
        balance=1-internal
        if j in deadA:
            A.append(0);B.append(balance)
        elif j in deadB:
            A.append(balance);B.append(0)
        else:
            coefficient=next(choices)
            A.append(coefficient);B.append(balance-coefficient)
    raise AssertionError('previous row must have both terms zero')


def basis(poly):
    return [int(poly.nth(i)) for i in range(6)]


def run():
    t,u=symbols('t u')
    H=Poly(t**6-4*t**5+6*t**4-7*t**3+4*t*t-2*t+1,t,domain=ZZ)
    left,right=Fraction(71942965,10**8),Fraction(71942966,10**8)
    assert H.count_roots(0,1)==1 and H.count_roots(Rational(left),Rational(right))==1
    binary=Poly(H.as_expr(),t,modulus=2)
    assert binary.is_irreducible
    for k,z in FORCED.items():
        assert implied_zeros(k)==z
    result={}
    for k,lo,hi in ((12,4877,4878),(13,17567,17569)):
        value=forced_previous(k,Interval(left,right),[])
        assert Fraction(lo,10**4)<value.lo<=value.hi<Fraction(hi,10**4)
        assert not value.lo<=1<=value.hi
        result[str(k)]={'forced_zeros':FORCED[k],'exact_interval':[str(value.lo),str(value.hi)],
                        'convenient_bounds':[str(Fraction(lo,10**4)),str(Fraction(hi,10**4))]}
    R=Poly(forced_previous(14,t,[u]),u)
    assert R.degree()==1
    base=Poly(R.nth(0),t,domain=ZZ).rem(H)
    slope=Poly(R.nth(1),t,domain=ZZ).rem(H)
    numerator=int(H.resultant(Poly(1,t,domain=ZZ)-base))
    denominator=int(H.resultant(slope))
    assert numerator==-8060939 and denominator==-208477
    assert numerator%denominator!=0
    result['14']={'forced_zeros':FORCED[14],'base':basis(base),'slope':basis(slope),
                  'numerator_norm':numerator,'denominator_norm':denominator,
                  'nonintegral_norm':str(Fraction(numerator,denominator))}
    w=symbols('w')
    low_A,low_B=lower(t)
    q=1-low_A[3]*low_B[5]
    affine_cases=((15,'free',[u,w],(31,1272)),
                  (16,'A8=0',[0,u,w],(1,261)),
                  (16,'A8=q',[q,u,w],(88,27)))
    for k,tag,splits,expected in affine_cases:
        residual=Poly(forced_previous(k,t,splits)-1,u,w)
        assert residual.total_degree()==1
        coefficients=[Poly(residual.coeff_monomial(g),t,domain=ZZ).rem(H) for g in (u,w,1)]
        assert basis(coefficients[0])==[-6,10,-16,21,-12,4]
        assert basis(coefficients[1])==[2,-8,21,-24,20,-6]
        checks=[]
        for (prime,root),constant in zip(((163,-6),(1279,168)),expected,strict=True):
            assert H.eval(root)%prime==0
            values=[int(f.eval(root))%prime for f in coefficients]
            assert values==[0,0,constant] and constant!=0
            checks.append({'prime':prime,'root':root,'values':values})
        result.setdefault(str(k),{'forced_zeros':FORCED[k],'cases':[]})['cases'].append(
            {'case':tag,'coefficients':[basis(f) for f in coefficients],'modular_checks':checks})
    return {'status':'verified','rows':result,'root_interval':[str(left),str(right)],
            'python':sys.version,'sympy':sympy_version}


def main():
    if not __debug__:
        raise SystemExit('assertions required')
    p=ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    result=run()
    with Path(__file__).open('rb') as f:
        result['source_sha256']=file_digest(f,'sha256').hexdigest()
    with a.out.open('x') as f:
        f.write(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print({k:v.get('convenient_bounds',v.get('nonintegral_norm')) for k,v in result['rows'].items()})


if __name__=='__main__':
    main()
