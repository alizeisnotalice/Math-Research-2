#!/usr/bin/env python3
"""Two short support exclusions, reconstructed using integer polynomials."""
import json
from verify_next_strip import add as uadd,mul as umul,scale as uscale


def add(*polys):
    out={}
    for p in polys:
        for e,v in p.items():out[e]=out.get(e,0)+v
    return {e:v for e,v in out.items() if v}


def scale(p,c):return {e:v*c for e,v in p.items() if v*c}


def mul(p,q):
    out={}
    for e,v in p.items():
        for f,w in q.items():
            g=tuple(a+b for a,b in zip(e,f));out[g]=out.get(g,0)+v*w
    return {e:v for e,v in out.items() if v}


one={(0,0,0):1};p={(1,0,0):1};q={(0,1,0):1};r={(0,0,1):1}


def prefix(I,J,A,last):
    support={i+j for i in I for j in J};B={}
    for k in range(last+1):
        B[k]=add(one if k in support else {},*(scale(mul(v,B.get(k-i,{})),-1)
                     for i,v in A.items() if i))
    return B


def qcoeff(poly,j):
    assert all(e[2]==0 for e in poly)
    selected={e[0]:v for e,v in poly.items() if e[1]==j}
    return [selected.get(i,0) for i in range(max(selected,default=0)+1)]


def check():
    I=[0,3,6,7,11];J=[0,3,6,7,9,10,14,15,19]
    B=prefix(I,J,{0:one,3:p,6:q,7:r,11:one},12)
    f=B[12]
    assert all(e[1]<=2 for e in f) and qcoeff(f,2)==[1]
    discriminant=uadd(umul(qcoeff(f,1),qcoeff(f,1)),uscale(qcoeff(f,0),-4))
    assert discriminant==umul([-1,1],[3,3,-3,5])
    assert 3 in I and 3 in J and 12 not in J
    I2=[0,1,2,4,10];J2=[0,2,3,4,5,7,11,13,19]
    # B1=0, c1=1 force A1=1 before this recurrence.
    B=prefix(I2,J2,{0:one,1:one,2:p,4:q,10:one},9)
    assert not B[1] and all(j not in J2 for j in (6,8,9))
    f6=B[6];assert all(e[1]<=1 for e in f6)
    a,b=qcoeff(f6,1),qcoeff(f6,0)
    def resultant(g):
        assert all(e[1]<=2 for e in g)
        return uadd(uadd(umul(umul(a,a),qcoeff(g,0)),
                         uscale(umul(umul(a,b),qcoeff(g,1)),-1)),
                    umul(umul(b,b),qcoeff(g,2)))
    P,R=resultant(B[8]),resultant(B[9])
    assert P==[-1,10,-23,24,-10,2,-1]
    assert R==[4,-20,42,-54,40,-14]
    U=[1701212,-4634648,7348502,-6423148,2455950]
    V=[526437,-2779507,3780591,-1638862,308432,-175425]
    assert uadd(umul(U,P),umul(V,R))==[404536]
    return {'verified':True,'arithmetic':'integer polynomials only','cases':[
        {'I':I,'J':J,'zero_index':12,'certificate':'monic quadratic in A6 has discriminant (p-1)(3+3p(1-p)+5p^3)<0 for 0<p=A3<1'},
        {'I':I2,'J':J2,'forced_A1':1,'zero_indices':[6,8,9],
         'P_ascending':P,'R_ascending':R,'U_ascending':U,'V_ascending':V,
         'bezout_constant':404536,'certificate':'U*Res_q(f6,f8)+V*Res_q(f6,f9)=404536'}]}


if __name__=='__main__':print(json.dumps(check(),indent=2))
