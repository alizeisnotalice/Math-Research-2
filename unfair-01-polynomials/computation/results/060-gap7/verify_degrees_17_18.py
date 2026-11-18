#!/usr/bin/env python3
"""Exact Sturm, rational bounds, and mod-163 saturation replay."""
from argparse import ArgumentParser
from dataclasses import dataclass
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from pathlib import Path
import json
import sys
import time

H=(1,-2,4,-7,6,-4,1)
P=163


def evaluate(poly,t):
    result=0
    for a in reversed(poly):
        result=result*t+a
    return result


def negative_remainder(a,b):
    a=list(a)
    while a and len(a)>=len(b):
        shift=len(a)-len(b)
        scale=a[-1]/b[-1]
        for j,v in enumerate(b):
            a[j+shift]-=scale*v
        while a and a[-1]==0:
            a.pop()
    return tuple(-z for z in a)


def sign(x):
    return (x>0)-(x<0)


def sign_changes(signs):
    p=[x for x in signs if x]
    return sum(x!=y for x,y in zip(p,p[1:]))


@dataclass(frozen=True)
class Interval:
    lo: F
    hi: F

    def __add__(self,other):
        other=to_interval(other)
        return Interval(self.lo+other.lo,self.hi+other.hi)

    __radd__=__add__

    def __neg__(self):
        return Interval(-self.hi,-self.lo)

    def __sub__(self,other):
        return self+-to_interval(other)

    def __rsub__(self,other):
        return to_interval(other)+-self

    def __mul__(self,other):
        other=to_interval(other)
        endpoints=[x*y for x,y in product((self.lo,self.hi),(other.lo,other.hi))]
        return Interval(min(endpoints),max(endpoints))

    __rmul__=__mul__


def to_interval(x):
    return x if isinstance(x,Interval) else Interval(F(x),F(x))


def common_values(t):
    m=1-t+t*t
    n=1-t*m
    p=1-(1-t)*n
    v=t*t*m
    q=1-n*v
    x=1-p*v-(1-t)*q
    y=1-p*v-t*q
    return m,n,p,v,q,x,y


def real_bounds():
    h=tuple(map(F,H))
    h1=tuple(j*h[j] for j in range(1,len(h)))
    sturm=[h,h1]
    while True:
        r=negative_remainder(*sturm[-2:])
        if not r:
            break
        sturm.append(r)
    signs=[[sign(evaluate(g,t)) for g in sturm] for t in (0,1)]
    assert list(map(sign_changes,signs))==[4,3]
    lo,hi=F(71942965,10**8),F(71942966,10**8)
    assert evaluate(H,lo)>0>evaluate(H,hi)
    t=Interval(lo,hi)
    m,n,p,v,q,x,y=common_values(t)
    delta=1-m*q-t*x
    assert lo>F(1,2)
    assert all(0<z.lo<=z.hi<1 for z in (m,n,p,v,q,x,y))
    assert delta.lo>F(508,10000)
    return {'sturm_signs_0_1':signs,'t_isolation':[str(lo),str(hi)],
            'strict_positive_y_and_delta_lower_bounds':[str(y.lo),str(delta.lo)]}


def affine_system(j,left8,left9):
    # Coordinates are the constant, the C11..C(j-1) bits, and the two free
    # tail coefficients. All lower-row products must have a constant factor.
    names=['1',*[f'C{k}' for k in range(11,j)],'u','w']
    count=len(names)
    zero=(0,)*count

    def const(x):
        return (x%P,)+zero[1:]

    def basis(name):
        k=names.index(name)
        return zero[:k]+(1,)+zero[k+1:]

    def add(u,v):
        return tuple((a+b)%P for a,b in zip(u,v,strict=True))

    def subtract(u,v):
        return tuple((a-b)%P for a,b in zip(u,v,strict=True))

    def left_times(n,u):
        return tuple(n*a%P for a in u)

    def multiply(u,v):
        if not any(u[1:]):
            return left_times(u[0],v)
        if not any(v[1:]):
            return left_times(v[0],u)
        raise AssertionError('unexpected nonlinear term')

    t=-6
    m,n,p,v,q,x,y=common_values(t)
    A=list(map(const,(1,t,0,n,p,0,0,0,q if left8 else 0,(x if left8 else y) if left9 else 0)))
    B=list(map(const,(1,1-t,m,0,0,v,0,0,0 if left8 else q,0 if left9 else (x if left8 else y))))
    free=iter(('u','w'))
    lockA={j-1,j-2,j-5}
    lockB={j-1,j-3,j-4}
    if j==18:
        (lockB if left8 else lockA).add(10)
    final=None
    for k in range(10,j):
        X=zero
        for i in range(1,k):
            X=add(X,multiply(A[i],B[k-i]))
        C=const(1) if k==10 else basis(f'C{k}')
        S=subtract(C,X)
        if k in lockA and k in lockB:
            assert k==j-1
            final=S
            A.append(zero)
            B.append(zero)
        elif k in lockA:
            A.append(zero)
            B.append(S)
        elif k in lockB:
            A.append(S)
            B.append(zero)
        else:
            a=basis(next(free))
            A.append(a)
            B.append(subtract(S,a))
    assert list(free)==[] and final is not None
    assert final[-2:]==(0,0)
    return names, final


def residue_filter():
    assert P==163 and all(P%i for i in range(2,13))
    assert evaluate(H,-6)%P==0
    common={
        17:(112,161,11,78,6,1),
        18:(120,112,161,11,78,6,1),
    }
    expected_constants={
        (17,False,False):14,(17,True,True):156,
        (18,False,False):46,(18,False,True):61,
        (18,True,False):47,(18,True,True):35,
    }
    total=0
    survivors=0
    table=[]
    for (j,left8,left9),expected in expected_constants.items():
        names,res=affine_system(j,left8,left9)
        constant=res[0]
        weights=res[1:-2]
        assert (constant,weights)==(expected,common[j])
        survivors_here=[]
        for digits in product((0,1),repeat=len(weights)):
            total+=1
            if (constant+sum(x*y for x,y in zip(weights,digits,strict=True)))%P==0:
                assert digits[0]==0
                survivors_here.append(digits)
        survivors+=len(survivors_here)
        table.append({'degree':j,'row8_side':'A' if left8 else 'B','row9_side':'A' if left9 else 'B',
                      'coordinates':names,'residue':res,'bit_patterns_remaining':survivors_here,
                      'compatible_with_C11_equal_1':0})
    assert (total,survivors)==(640,4)
    return {'modulus':P,'t_residue':-6,'total_binary_patterns':total,
            'residue_patterns_before_C11':survivors,'live':0,'cases':table}


def main():
    if not __debug__:raise SystemExit('run with assertions enabled')
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    arg=parser.parse_args()
    start=time.monotonic()
    data={'status':'verified','real_bounds':real_bounds(),'residues':residue_filter(),
          'python':sys.version,'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
          'seconds':time.monotonic()-start}
    with arg.out.open('x') as f:f.write(json.dumps(data,sort_keys=True,indent=2)+'\n')
    print(json.dumps({'status':data['status'],'residues':{k:v for k,v in data['residues'].items() if k!='cases'}},sort_keys=True))


if __name__=='__main__':
    main()
