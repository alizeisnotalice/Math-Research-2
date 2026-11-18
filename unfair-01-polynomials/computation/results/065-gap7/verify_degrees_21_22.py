"""Exact saturation certificate for factor degrees 21 and 22."""
from argparse import ArgumentParser
from collections import Counter
from hashlib import sha256
from itertools import product, permutations
from fractions import Fraction as Q
from pathlib import Path
import json
import sys
import time

sys.dont_write_bytecode=True
from verify_degree_19 import Element, norm, cross, element
from verify_degrees_17_18 import P, H, Interval, real_bounds, evaluate, common_values

T=Element((0,1,0,0,0,0))
ZERO=element(0)


def setup(t):
    m,n,p,v,_,_,_=common_values(t)
    return [1,t,0,n,p,0,0,0],[1,1-t,m,0,0,v,0,0]


def split(s,typ,h=0):
    if typ=='A':return s,0
    if typ=='B':return 0,s
    if typ=='X':return h,s-h
    assert typ=='0'
    return 0,0


def lock(j,early,middle,bits):
    sequence=list(early)
    if j%2==0:sequence.append(middle)
    sequence.extend(('0' if p=='X' else p) for p in early[::-1])
    assert 8+len(sequence)-1==j-8
    cs=[1,1,1,1]+list(bits)
    assert len(sequence)==len(cs)
    return sequence,cs


def frame(t,j,early,middle,bits,free):
    A,B=setup(t)
    residuals=[]
    supported=iter(free)
    seq,codes=lock(j,early,middle,bits)
    for k,(kind,C) in enumerate(zip(seq,codes,strict=True),start=8):
        s=C-cross(A,B,k)
        var=next(supported) if kind=='X' else 0
        a,b=split(s,kind,var)
        A.append(a);B.append(b)
        if kind=='0':residuals.append(s)
    assert not list(supported)
    return residuals,A,B


def positive_low_rows():
    return ({0,1,3,4},{0,1,2,5})


def impossible_bits(j,early,bits):
    ap,bp=positive_low_rows()
    for i,position in enumerate(early,start=8):
        if position in ('A','X'):ap.add(i)
        if position in ('B','X'):bp.add(i)
    return any(not c and any(k-a in bp for a in ap if 1<=a<k) for k,c in enumerate(bits,start=12))


def determinant(A):
    if any(len(r)!=len(A) for r in A):raise ValueError('square matrix required')
    s=0
    for p in permutations(range(len(A))):
        sign=(-1)**sum(x>y for i,x in enumerate(p) for y in p[i+1:])
        value=sign
        for i,j in enumerate(p):
            value=value*A[i][j]
        s=s+value
    return s


def solve(t,j,early,middle,bits):
    dimension=early.count('X')
    zero=(0,)*dimension
    base=frame(t,j,early,middle,bits,zero)[0]
    assert len(base)==dimension
    rhs=[-v for v in base]
    cols=[]
    for k in range(dimension):
        x=zero[:k]+(1,)+zero[k+1:]
        cols.append([a-b for a,b in zip(frame(t,j,early,middle,bits,x)[0],base,strict=True)])
    matrix=list(zip(*cols,strict=True))
    det=determinant(matrix)
    if det==ZERO:return ('singular',matrix,rhs)
    nd=norm(det)
    assert nd
    for i in range(dimension):
        other=matrix[:]
        other=[r[:i]+(b,)+r[i+1:] for r,b in zip(matrix,rhs,strict=True)]
        nn=norm(determinant(other))
        ratio=Q(nn,nd)
        if ratio.denominator>1:return ('norm',i,nn,nd)
    return ('uncertain_norm',matrix,rhs)


def signs_for_pure(j,early,middle,bits):
    lo,hi=map(Q,real_bounds()['t_isolation'])
    t=Interval(lo,hi)
    _,A,B=frame(t,j,early,middle,bits,[])
    for k in range(8,len(A)):
        for side,arr in (('A',A),('B',B)):
            v=arr[k]
            if isinstance(v,Interval) and v.hi<0:
                return {'nonnegative_rejection':(side,k),'rational_upper_bound':str(v.hi)}
    return None


def proof_exists_positive(A,B,k):
    return any(a+b==k for a in A for b in B)


def precise_support(j,early,middle,bits):
    lo,hi=map(Q,real_bounds()['t_isolation'])
    t=Interval(lo,hi)
    _,A,B=frame(t,j,early,middle,bits,[])
    ap,bp=positive_low_rows()
    ap.update(i for i,v in enumerate(A) if isinstance(v,Interval) and v.lo>0)
    bp.update(i for i,v in enumerate(B) if isinstance(v,Interval) and v.lo>0)
    return ap,bp


def exclude_early(j):
    counters=Counter()
    undecided=[]
    pure=[]
    details=[]
    for early_seq in product('ABX',repeat=3):
        early=''.join(early_seq)
        for middle in ('AB' if j%2==0 else '') or (None,):
            for bits in product((0,1),repeat=j-19):
                counters['patterns']+=1
                tag=(early,middle,bits)
                if impossible_bits(j,early,bits):
                    counters['positive_forces_digit']+=1
                    continue
                if 'X' in early:
                    verdict,*data=solve(T,j,early,middle,bits)
                    counters[verdict]+=1
                    if verdict!='norm':undecided.append((tag,verdict))
                    else:details.append({'shape':tag,'norm_witness':data})
                else:
                    bound=signs_for_pure(j,early,middle,bits)
                    if bound:
                        counters['pure_negative']+=1
                        details.append({'shape':tag,**bound})
                    else:
                        counters['pure_possible']+=1
                        pure.append(tag)
    return {'j':j,'count':dict(counters),'undecided_mixed':undecided,'pure':pure,'details':details}


def residue(j,early,middle,bits,digits,free_values):
    _,A,B=frame(-6,j,early,middle,bits,[])
    A=[v%P for v in A]
    B=[v%P for v in B]
    lock_a={j-1,j-2,j-5}
    lock_b={j-1,j-3,j-4}
    values=iter(free_values)
    for k,digit in enumerate(digits,start=j-7):
        s=(digit-cross(A,B,k))%P
        if k in lock_a and k in lock_b:
            assert k==j-1 and not list(values)
            return s
        if k in lock_a:
            a,b=0,s
        elif k in lock_b:
            a,b=s,0
        else:
            a=next(values)
            b=(s-a)%P
        A.append(a)
        B.append(b)
    raise ValueError('missing last residual')


def exclude_late(j,prefixes):
    seqs=list(product((0,1),repeat=7))
    weight=(120,112,161,11,78,6,1)
    expected={
        (21,'AAB',None):(135,(14,15,17,18,19,20),((0,1,0,0,1,0,1),)),
        (21,'ABB',None):(101,(14,15,16,17,18,19,20),()),
        (22,'AAA','B'):(42,(15,16,17,18,19,20,21),((0,1,1,1,0,0,0),(1,0,0,0,0,0,1))),
        (22,'AAB','B'):(16,(15,16,18,19,20,21),((1,1,0,0,1,0,0),)),
        (22,'ABA','B'):(154,(15,16,17,19,21),((0,0,1,1,0,0,0),)),
        (22,'ABB','B'):(148,(15,16,17,18,19,20,21),((0,0,1,1,0,1,0),)),
    }
    assert (P,evaluate(H,-6)%P)==(163,0) and all(P%i for i in range(2,13))
    counter=0
    survivors=[]
    for early,middle,bits in prefixes:
        base=residue(j,early,middle,bits,seqs[0],(0,0))
        ap,bp=precise_support(j,early,middle,bits)
        forced=tuple(k for k in range(j-7,j) if proof_exists_positive(ap,bp,k))
        certain=[]
        live=[]
        for digits in seqs:
            counter+=1
            reduced=(base+sum(a*b for a,b in zip(weight,digits,strict=True)))%P
            for uv in ((0,0),(1,0),(0,1)):
                assert reduced==residue(j,early,middle,bits,digits,uv)
            if not reduced:
                certain.append(digits)
                if all(digits[k-(j-7)] for k in forced):
                    live.append(digits)
        assert (base,forced,tuple(certain))==expected[j,early,middle] and not live
        survivors.append({'shape':(early,middle,bits),'constant_mod_163':base,
                          'forced_ones':forced,'binary_patterns_with_zero_residue':certain})
    return {'tested_late_codes':counter,'common_weights':weight,'pure_shapes':survivors}


def main():
    if not __debug__:raise SystemExit('run with assertions enabled')
    ap=ArgumentParser(description=__doc__)
    ap.add_argument('--out',type=Path,required=True)
    opts=ap.parse_args()
    start=time.monotonic()
    observed=[]
    expected={
        21:dict(patterns=108,positive_forces_digit=79,pure_negative=8,pure_possible=2,norm=19),
        22:dict(patterns=432,positive_forces_digit=368,pure_negative=20,pure_possible=4,norm=40),
    }
    for j in (21,22):
        r=exclude_early(j)
        assert r['count']==expected[j] and not r['undecided_mixed']
        r['late']=exclude_late(j,r.pop('pure'))
        assert r['late']['tested_late_codes']==128*r['count']['pure_possible']
        observed.append(r)
    names=(Path(__file__).name,'verify_degree_19.py','verify_degrees_17_18.py')
    base=Path(__file__).parent
    report={'status':'verified','bounds':real_bounds(),'degrees':observed,
            'python':sys.version,'seconds':time.monotonic()-start,
            'source_sha256':{name:sha256((base/name).read_bytes()).hexdigest() for name in names}}
    with opts.out.open('x') as f:f.write(json.dumps(report,sort_keys=True,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'shape_counts':{r['j']:r['count'] for r in observed},
                      'late_modular_codes':sum(r['late']['tested_late_codes'] for r in observed)},sort_keys=True))


if __name__=='__main__':main()
