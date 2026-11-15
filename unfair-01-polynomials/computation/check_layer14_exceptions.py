#!/usr/bin/env python3
"""Independent exact checks for two non-domination cases in round 14."""
import json
from verify_next_strip import add,mul,scale


def check():
    I=[0,2,4,7]
    J=[0,2,6,8,9,11,12,14,17]
    sums={i+j for i in I for j in J}
    # c2=a+B2=1 forces 0<a<1. Absent B4 solves A4=a^2-a+1.
    A={0:[1],2:[0,1],4:[1,-1,1],7:[1]}
    B={}
    for k in range(11):
        p=[int(k in sums)]
        for i,v in A.items():
            if i:p=add(p,scale(mul(v,B.get(k-i,[0])),-1))
        B[k]=p
    assert B[4]==[0]
    assert 4 not in J and 10 not in J and 2 in J
    assert B[10]==mul([1,-2,1],[1,-1,1])
    # The other case is a certified specialization of the round-11 theorem.
    I2=[0,1,3,8]
    J2=[0,1,2,3,4,6,9,11,16]
    assert I2==[0,1,3,2*3+2] and {1,3}<=set(J2)
    return {'verified':True,'cases':[
        {'I':I,'J':J,'certificate':'B10=(a-1)^2*(a^2-a+1)>0, but 10 not in J',
         'domain':'0<a<1 from c2=a+B2=1'},
        {'I':I2,'J':J2,'certificate':'round 11, r=3, d=1, m=8',
         'domain':'0<A1<1 from c1; 0<A3<1 from c3=A3+A1*B2+B3 and B3>0',
         'dependency':'proofs/011-next-four-term-strip.md'}]}


if __name__=='__main__':print(json.dumps(check(),indent=2))
