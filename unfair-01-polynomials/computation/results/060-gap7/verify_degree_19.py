#!/usr/bin/env python3
"""Exact unbounded-extension exclusion of either factor degree 19."""
from argparse import ArgumentParser
from dataclasses import dataclass
from hashlib import sha256
from itertools import product
from pathlib import Path
import json
import sys
import time

sys.dont_write_bytecode=True
from verify_degrees_17_18 import H, P, F, Interval, common_values, evaluate, real_bounds


@dataclass(frozen=True)
class Element:
    """Integer element of the basis 1,t,...,t^5 in Z[t]/H."""
    coefficients: tuple[int,...]

    def __add__(self,other):
        other=element(other)
        return Element(tuple(a+b for a,b in zip(self.coefficients,other.coefficients,strict=True)))

    __radd__=__add__

    def __neg__(self):
        return Element(tuple(-a for a in self.coefficients))

    def __sub__(self,other):
        return self+-element(other)

    def __rsub__(self,other):
        return element(other)+-self

    def __mul__(self,other):
        other=element(other)
        out=[0]*11
        for j,a in enumerate(self.coefficients):
            for k,b in enumerate(other.coefficients):
                out[j+k]+=a*b
        for j in range(10,5,-1):
            for i,h in enumerate(H[:-1]):
                out[j-6+i]-=out[j]*h
        return Element(tuple(out[:6]))

    __rmul__=__mul__


def element(x):
    if isinstance(x,Element):return x
    return Element((x,0,0,0,0,0))


def determinant(table):
    a=[list(map(F,row)) for row in table]
    out=1
    for k in range(len(a)):
        j=next((i for i in range(k,len(a)) if a[i][k]),None)
        if j is None:return 0
        if j!=k:
            a[k],a[j]=a[j],a[k]
            out=-out
        lead=a[k][k]
        out*=lead
        for i in range(k+1,len(a)):
            h=a[i][k]/lead
            for j in range(k+1,len(a)):
                a[i][j]-=h*a[k][j]
    assert out.denominator==1
    return int(out)


def norm(x):
    unit=tuple(1 if j==k else 0 for j in range(6) for k in range(6))
    columns=[(x*Element(unit[k*6:(k+1)*6])).coefficients for k in range(6)]
    return determinant(zip(*columns,strict=True))


def prefix(t,A8,left9):
    m,n,p,v,q,x,y=common_values(t)
    B8=q-A8
    s9=1-p*v-t*B8-(1-t)*A8
    A9,B9=(s9,0) if left9 else (0,s9)
    s10=1-m*A8-t*B9-(1-t)*A9
    A10,B10=(s10,0) if left9 else (0,s10)
    A=[1,t,0,n,p,0,0,0,A8,A9,A10]
    B=[1,1-t,m,0,0,v,0,0,B8,B9,B10]
    return A,B


def cross(A,B,k):
    return sum(A[i]*B[k-i] for i in range(1,k))


def norm_rejections():
    t=Element((0,1,0,0,0,0))
    expected={True:(-359,1171),False:(587,-1)}
    result={}
    for left in (False,True):
        A,B=prefix(t,0,left)
        const=cross(A,B,11)-1
        A,B=prefix(t,1,left)
        slope=cross(A,B,11)-1-const
        denominator,numerator=norm(slope),norm(-const)
        assert (denominator,numerator)==expected[left]
        value=F(numerator,denominator)
        assert value.denominator>1
        result['A' if left else 'B']={'alpha':slope.coefficients,'gamma':(-const).coefficients,
                                      'norm_alpha':denominator,'norm_gamma':numerator,'norm_A8':str(value)}
    return result


def remainder(left8,left9,digits,free):
    t=-6
    q=common_values(t)[4]
    A,B=prefix(t,q if left8 else 0,left9)
    A=[a%P for a in A];B=[b%P for b in B]
    s=(1-cross(A,B,11))%P
    A.append(s if left8 else 0)
    B.append(0 if left8 else s)
    lockA={18,17,14}
    lockB={18,16,15}
    free=iter(free)
    for k,c in enumerate(digits,start=12):
        total=(c-cross(A,B,k))%P
        if k in lockA and k in lockB:
            assert k==18 and list(free)==[]
            return total
        if k in lockA:
            A.append(0);B.append(total)
        elif k in lockB:
            A.append(total);B.append(0)
        else:
            a=next(free)
            A.append(a);B.append((total-a)%P)
    raise ValueError('missing row-18 residual')


def verify():
    bounds=real_bounds()
    lo,hi=map(F,bounds['t_isolation'])
    t=Interval(lo,hi)
    q=common_values(t)[4]
    assert P==163 and all(P%i for i in range(2,13)) and evaluate(H,-6)%P==0
    expected_constants={(False,False):117,(False,True):162,(True,False):121,(True,True):131}
    weights=(120,112,161,11,78,6,1)
    exact=[]
    all_digits=list(product((0,1),repeat=7))
    for pair,constant in expected_constants.items():
        a,b=pair
        A,B=prefix(t,q if a else 0,b)
        c=cross(A,B,11)
        if pair==(False,False):
            assert c.lo>1
        else:
            assert 0<c.lo<c.hi<1
        zero=all_digits[0]
        assert all(remainder(a,b,zero,f)==constant for f in ((0,0),(1,0),(0,1)))
        survivors=[]
        for digits in all_digits:
            value=(constant+sum(x*y for x,y in zip(digits,weights,strict=True)))%P
            assert value==remainder(a,b,digits,(0,0))
            if not value:survivors.append(digits)
        eligible=[]
        for d in survivors:
            if not d[0]:
                continue
            if pair==(False,False):
                continue
            if pair==(False,True) and d[2]==0:
                continue
            if pair==(True,False) and d[1]==0:
                continue
            eligible.append(d)
        assert not eligible
        exact.append({'row8':'A' if a else 'B','row9':'A' if b else 'B',
                      'X11_rational_interval':[str(c.lo),str(c.hi)],
                      'residual_constant':constant,'weights_C12_to_C18':weights,
                      'residue_survivors':survivors,'eligible':eligible})
    return {'result':'verified','real_bounds':bounds,'modular_rows':exact,'modulus':P,'t_residue':-6,
            'total_bit_patterns':4*len(all_digits),'A8_norm_rejections':norm_rejections()}


def main():
    if not __debug__:raise SystemExit('run with assertions enabled')
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--out',required=True,type=Path)
    args=parser.parse_args()
    start=time.monotonic()
    result=verify()
    sources=[Path(__file__),Path(__file__).with_name('verify_degrees_17_18.py')]
    result.update(python=sys.version,seconds=time.monotonic()-start,
                  sources_sha256={p.name:sha256(p.read_bytes()).hexdigest() for p in sources})
    with args.out.open('x') as out:out.write(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('result','total_bit_patterns','A8_norm_rejections')},sort_keys=True))


if __name__=='__main__':
    main()
