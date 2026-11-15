#!/usr/bin/env python3
"""Independent integer-polynomial certificate; no CAS or searcher imported."""
import json


def add(*polys):
    r = {}
    for p in polys:
        for e,v in p.items():
            r[e] = r.get(e,0)+v
    return {e:v for e,v in r.items() if v}


def scale(p,n):
    return {e:n*v for e,v in p.items() if n*v}


def mul(p,q):
    r = {}
    for (i,j),v in p.items():
        for (k,l),w in q.items():
            e = (i+k,j+l)
            r[e] = r.get(e,0)+v*w
    return {e:v for e,v in r.items() if v}


def check():
    I = [0,1,4,5,11]
    J = [0,1,2,3,5,8,9,15]
    sums = {i+j for i in I for j in J}
    assert set(range(11)) <= sums
    one,a,c = {(0,0):1},{(1,0):1},{(0,1):1}
    # b=a^4-a^3+a^2-a+1 follows from the absent B coefficient at 4.
    b = {(i,0):(-1)**i for i in range(5)}
    A = {0:one,1:a,4:b,5:c,11:one}
    B = {}
    for k in range(8):
        B[k] = add(one, *(scale(mul(v,B.get(k-i,{})),-1)
                           for i,v in A.items() if i))
    assert not B[4]
    g6,g7 = B[6],B[7]
    assert all(j <= 1 for _,j in g6.keys() | g7.keys())
    l6 = {(i,0):v for (i,j),v in g6.items() if j == 1}
    l7 = {(i,0):v for (i,j),v in g7.items() if j == 1}
    assert l6 == {(1,0):2,(0,0):-1}
    assert l7 == {(2,0):-3,(1,0):2,(0,0):-1}
    resultant = add(mul(l6,g7),scale(mul(l7,g6),-1))
    expected = mul(mul(mul({(3,0):1},add(a,scale(one,-2))),
                       add(a,scale(one,-1))),add(mul(a,a),one))
    assert resultant == expected
    return {'verified':True,'I':I,'J':J,'zero_B_indices':[4,6,7],
            'identity':'(2a-1)g7-(-3a^2+2a-1)g6=a^3(a-2)(a-1)(a^2+1)',
            'strict_domain':'0<a<1, since a+B1=1 and both supports contain 1',
            'contradiction':'left=0; right>0',
            'coefficient_arithmetic':'Python integers only'}


if __name__ == '__main__':
    print(json.dumps(check(),indent=2))
