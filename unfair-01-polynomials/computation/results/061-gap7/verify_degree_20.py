#!/usr/bin/env python3
"""Exact certification of all saturation shapes for candidate factor degree 20."""
from argparse import ArgumentParser
from hashlib import sha256
from itertools import product
from pathlib import Path
import json
import sys
import time

sys.dont_write_bytecode=True
from verify_degree_19 import Element, cross, norm
from verify_degrees_17_18 import H, P, F, Interval, common_values, evaluate, real_bounds


def placement(value,side,w=0):
    if side=='A':return value,0
    if side=='B':return 0,value
    assert side=='X'
    return w,value-w


def frame(t,shape,u=0,w=0):
    s8,s9,s10=shape
    m,n,p,v,q,x,y=common_values(t)
    A=[1,t,0,n,p,0,0,0]
    B=[1,1-t,m,0,0,v,0,0]
    a,b=placement(q,s8,u)
    A.append(a);B.append(b)
    a,b=placement(1-cross(A,B,9),s9,w)
    A.append(a);B.append(b)
    a,b=placement(1-cross(A,B,10),s10)
    A.append(a);B.append(b)
    x11=cross(A,B,11)
    a,b=(0,0) if s9=='X' else placement(1-x11,s9)
    A.append(a);B.append(b)
    return A,B,x11,cross(A,B,12)


def mixed_residual(t,shape,uv):
    _,_,x11,x12=frame(t,shape,*uv)
    out=[]
    if shape[1]=='X':out.append(x11-1)
    if shape[0]=='X':out.append(x12-1)
    return out


def mixed():
    t=Element((0,1,0,0,0,0))
    expected={
        'AXA':(5855,11),'AXB':(21289,379),'BXA':(-31139,11),'BXB':(-1,379),
        'XAA':(274959,32960),'XAB':(-1486327,20407),'XBA':(206015,74569),'XBB':(-17,31),
        'XXA':(4807601,333325),'XXB':(925,7603),
    }
    out={}
    for shape,wanted in expected.items():
        cols=['u'] if shape[0]=='X' else []
        if shape[1]=='X':cols.append('w')
        const=mixed_residual(t,shape,(0,0))
        lhs=[]
        for variable in cols:
            coord=(1,0) if variable=='u' else (0,1)
            lhs.append([a-b for a,b in zip(mixed_residual(t,shape,coord),const,strict=True)])
        rhs=[-r for r in const]
        if len(cols)==1:
            denominator=lhs[0][0]
            numerator=rhs[0]
        else:
            (a,c),(b,d)=lhs
            e,f=rhs
            denominator=a*d-b*c
            numerator=e*d-b*f
        nd,nn=norm(denominator),norm(numerator)
        assert nd
        ratio=F(nn,nd)
        assert ratio==F(*wanted) and ratio.denominator>1
        out[shape]={'unique_coefficient':cols[0],'norm_numerator':nn,'norm_denominator':nd,'norm_quotient':str(ratio)}
    assert len(out)==10
    return out


def residue(sides,digits,free):
    A,B,_,x12=frame(-6,sides)
    a,b=placement(1-x12,sides[0])
    A.append(a);B.append(b)
    A=[k%P for k in A];B=[k%P for k in B]
    free=iter(free)
    for k,d in enumerate(digits,start=13):
        residual=(d-cross(A,B,k))%P
        if k==19:
            assert not list(free)
            return residual
        if k in (18,15):
            a,b=0,residual
        elif k in (17,16):
            a,b=residual,0
        else:
            a=next(free)
            b=(residual-a)%P
        A.append(a);B.append(b)
    raise ValueError('missing residual')


def pure():
    data=real_bounds()
    lo,hi=map(F,data['t_isolation'])
    t=Interval(lo,hi)
    assert P==163 and all(P%i for i in range(2,13)) and evaluate(H,-6)%P==0
    coeffs=(120,112,161,11,78,6,1)
    constants={'AAA':109,'AAB':159,'ABA':110,'ABB':143}
    cases={}
    for sides in map(''.join,product('AB',repeat=3)):
        _,_,x11,x12=frame(t,sides)
        item={}
        if sides[0]=='B':
            if x11.lo>1:
                item['reject']='X11>1'
            else:
                assert 0<x11.lo<=x11.hi<1 and x12.lo>1
                item['reject']='X12>1'
        else:
            assert 0<x11.lo<=x11.hi<1 and 0<x12.lo<=x12.hi<1
            constant=constants[sides]
            zero=(0,)*7
            assert all(residue(sides,zero,f)==constant for f in ((0,0),(1,0),(0,1)))
            survivors=[];live=[]
            for digits in product((0,1),repeat=7):
                r=(constant+sum(a*b for a,b in zip(coeffs,digits,strict=True)))%P
                assert r==residue(sides,digits,(0,0))
                if r:continue
                survivors.append(digits)
                if not digits[0]:continue
                if sides[2]=='A' and not digits[2]:continue
                live.append(digits)
            assert not live
            item.update(residue_constant=constant,residue_C13_to_C19=coeffs,survivors=survivors,live=live)
        item['X11_interval']=[str(x11.lo),str(x11.hi)]
        item['X12_interval']=[str(x12.lo),str(x12.hi)]
        cases[sides]=item
    return data,cases


def main():
    if not __debug__:raise SystemExit('assertions required')
    ap=ArgumentParser(description=__doc__)
    ap.add_argument('--out',required=True,type=Path)
    args=ap.parse_args()
    start=time.monotonic()
    a=mixed()
    b,c=pure()
    names=('verify_degrees_17_18.py','verify_degree_19.py','verify_degree_20.py')
    source=Path(__file__).parent
    result={'status':'verified','real_bounds':b,'modulus':163,'t_residue':-6,
            'mixed_cases':a,'pure_cases':c,'total_saturation_shapes':18,'modular_binary_patterns':512,
            'python':sys.version,'seconds':time.monotonic()-start,
            'sources_sha256':{n:sha256((source/n).read_bytes()).hexdigest() for n in names}}
    with args.out.open('x') as f:f.write(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('status','total_saturation_shapes','modular_binary_patterns','seconds')},sort_keys=True))


if __name__=='__main__':
    main()
