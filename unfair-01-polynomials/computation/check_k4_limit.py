#!/usr/bin/env python3
"""Rebuild all eight k=4 limiting polynomials with integer arithmetic."""
from fractions import Fraction as Q
from itertools import product
import json
from verify_next_strip import add, mul, scale, sturm, count, sign_at, variations, W


def check():
    one, s, h = [1], [0,1], [1,1]
    w = W(4)
    h2w2 = mul(mul(h,h),mul(w,w))
    denominator = mul(h,h2w2)
    Z3 = [1,-3,6,-10]
    # E_all = 1/h - 1/(h^3 W4) + Z3/(h^2 W4^2).
    base = add(add(h2w2,scale(w,-1)),mul(h,Z3))
    assert base == [1,-2,3,-4,-10,6,-7,8,-9,10,25]
    lo, hi = Q(4,5),Q(6,7)
    iso = [Q(975853181,1190734317),Q(1136568134,1386838417)]
    # P2-1 = 2s^2-s + tau*(beta6-s*beta5).
    overflow = add(add(mul([0,-1,2],h2w2),mul([1,-1],w)),
                   add(mul(s,W(5)),scale(W(6),-1)))
    assert lo < iso[0] < iso[1] < hi
    rows = []
    for u,v,z in product((0,1),repeat=3):
        d0 = [u-1]
        d1 = add([v-1],scale(mul(s,d0),-1))
        d2 = add([z-1],scale(mul(s,d1),-1))
        p = add(base,scale(mul(mul(s,d2),denominator),-1))
        seq = sturm(p)
        roots = count(seq,lo,hi)
        assert roots == int((u,v,z)==(1,1,1))
        row = {'branch':[u,v,z],'polynomial_ascending':p,
               'endpoint_signs':[sign_at(p,lo),sign_at(p,hi)],
               'variations':[variations(seq,lo),variations(seq,hi)],'roots':roots}
        if roots:
            assert count(seq,*iso)==1
            qseq=sturm(overflow)
            assert count(qseq,*iso)==0
            assert sign_at(overflow,sum(iso)/2)==1
            row.update(isolator=list(map(str,iso)),overflow_numerator_ascending=overflow,
                       overflow_roots_in_isolator=0,overflow_sign=1)
        rows.append(row)
    return {'verified':True,'arithmetic':'integer primitive Sturm and Fraction endpoints',
            'range':list(map(str,(lo,hi))),'branches':rows,
            'conclusion':'All eight limiting branches excluded; no explicit R4.'}


if __name__=='__main__':
    print(json.dumps(check(),indent=2))
