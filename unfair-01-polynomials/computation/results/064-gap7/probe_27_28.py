"""Complete early prefix census with sparse polynomial arithmetic in Q(t)."""
from collections import Counter
from dataclasses import dataclass
from itertools import product
from pathlib import Path
import json,sys

sys.dont_write_bytecode=True
from verify_degrees_25_26 import (T,ZERO,INTERVAL,setup,positive_low_rows,
                                  proof_exists_positive,linear_test,local_linear_rejection)
from verify_degree_19 import cross,element

CODES='ABXZ'
PARTNER={'A':'A','B':'B','X':'Z','Z':'F'}


def initial_shapes(j):
    for states in product('ABX','ABX','ABX',CODES,CODES,CODES):
        early=''.join(states)
        for mid in ('AB' if j%2==0 else (None,)):
            yield early,mid


def required_positive(early):
    a,b=positive_low_rows()
    for k,state in enumerate(early,8):
        if state in ('A','X'):a.add(k)
        if state in ('B','X'):b.add(k)
    return a,b


def masks(j,early):
    a,b=required_positive(early)
    return [((1,) if proof_exists_positive(a,b,k) else (0,1)) for k in range(12,j-7)]


def sequence(j,shape):
    early,mid,_=shape
    a=early+(mid or '')+''.join(PARTNER[z] for z in early[::-1])
    assert len(a)==j-15
    return a


def slots(j,shape):
    return tuple(k for k,z in enumerate(sequence(j,shape),8) if z in ('X','F'))


def model(t,j,shape,free):
    values=iter(free)
    residual=[]
    A,B=setup(t)
    for k,(side,c) in enumerate(zip(sequence(j,shape),[1]*4+list(shape[2]),strict=True),8):
        x=c-cross(A,B,k)
        if side=='A':a,b=x,0
        elif side=='B':a,b=0,x
        elif side in ('X','F'):
            a=next(values);b=x-a
        else:
            assert side=='Z'
            a=b=0
            residual.append((k,x))
        A.append(a);B.append(b)
    assert not list(values)
    return residual,A,B


@dataclass(frozen=True)
class Polynomial:
    coef:dict

    def __add__(self,r):
        r=polynomial(r)
        z=self.coef.copy()
        for a,b in r.coef.items():
            z[a]=z.get(a,ZERO)+b
            if z[a]==ZERO:
                z.pop(a)
        return Polynomial(z)

    __radd__=__add__

    def __neg__(self):
        return Polynomial({a:-b for a,b in self.coef.items()})

    def __sub__(self,r):return self+-polynomial(r)

    def __rsub__(self,r):return polynomial(r)+-self

    def __mul__(self,r):
        r=polynomial(r)
        z={}
        for a,b in self.coef.items():
            for c,d in r.coef.items():
                k=tuple(sorted(a+c))
                assert len(k)<=2
                z[k]=z.get(k,ZERO)+b*d
        return Polynomial({a:b for a,b in z.items() if b!=ZERO})

    __rmul__=__mul__


def polynomial(x):
    if isinstance(x,Polynomial):return x
    e=element(x)
    return Polynomial({():e} if e!=ZERO else {})


def exact_system(j,shape):
    n=len(slots(j,shape))
    free=[Polynomial({(i,):element(1)}) for i in range(n)]
    rows,_,_=model(polynomial(T),j,shape,free)
    assert len(rows)==n
    return [(k,polynomial(v).coef) for k,v in rows]


def scout(j):
    count=Counter()
    late=[]
    early_cert=[]
    N=2**(j-19)
    for word,mid in initial_shapes(j):
        D=masks(j,word)
        count['all']+=N
        count['zeros']+=N-2**sum(len(d)==2 for d in D)
        for bits in product(*D):
            shape=(word,mid,bits)
            array=slots(j,shape)
            if not array:
                _,A,B=model(INTERVAL,j,shape,[])
                neg=next(((side,k,str(x.hi)) for side,row in (('A',A),('B',B)) for k,x in enumerate(row) if hasattr(x,'hi') and x.hi<0),None)
                if neg:
                    count['negative_pure']+=1
                    early_cert.append((shape,'negative_pure',neg))
                else:
                    count['pure_late']+=1
                    late.append(shape)
                continue
            E=exact_system(j,shape)
            linear=[f for _,f in E if all(len(m)<2 for m in f)]
            n=len(array)
            M=[tuple(f.get((i,),ZERO) for i in range(n)) for f in linear]
            rhs=[-f.get((),ZERO) for f in linear]
            if len(linear)==n:
                witness=linear_test(array,M,rhs)
                if witness:
                    count['norm']+=1
                    early_cert.append((shape,'norm',witness))
                    continue
            quick=local_linear_rejection(M,rhs)
            if quick:
                count[quick[0]]+=1
                early_cert.append((shape,*quick))
                continue
            count['systems_late']+=1
            late.append(shape)
    return {'j':j,'counts':dict(count),'late':late,'early_cert':early_cert}


if __name__=='__main__':
    if not __debug__:raise SystemExit('run with assertions')
    result=[]
    for j in map(int,sys.argv[1:]):
        a=scout(j)
        result.append(a)
        print(j,a['counts'],flush=True)
    Path('scout' + '.'.join(sys.argv[1:]) + '.json').write_text(json.dumps(result,separators=(',',':'))+'\n')
