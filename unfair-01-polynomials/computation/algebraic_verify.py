#!/usr/bin/env python3
"""Exact real-number-field factor verifier. Rational strings only, never floats.

Input: minpoly (ascending rational coefficients), isolating open interval,
A and B (dense coefficient lists; each coefficient is a polynomial in theta).
Irreducibility and root isolation are checked, not assumed. Sign evaluation uses
rational interval Horner arithmetic and Sturm root counts through SymPy.
"""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import sympy as sp


def rational(v):
    if not isinstance(v,(str,int)) or isinstance(v,bool):
        raise ValueError('Only exact rational strings or integers are accepted')
    return F(v)


class RealField:
    def __init__(self, polynomial, interval):
        self.x = sp.Symbol('t')
        self.p = sp.Poly(sum(sp.Rational(str(rational(c)))*self.x**i for i,c in enumerate(polynomial)), self.x, domain=sp.QQ)
        if self.p.degree() < 1 or not self.p.is_irreducible:
            raise ValueError('Minimal polynomial must be irreducible over Q')
        self.modulus = [F(str(c)) for c in reversed(self.p.monic().all_coeffs())]
        self.lo,self.hi = map(rational,interval)
        if not self.lo < self.hi or self.p.eval(self.lo)==0 or self.p.eval(self.hi)==0 or self.p.count_roots(self.lo,self.hi)!=1:
            raise ValueError('Interval must isolate one root, excluding endpoints')
        self.zero = self.norm([0])
        self.one = self.norm([1])

    def norm(self, coefficients):
        a = list(map(rational,coefficients))
        d = len(self.modulus)-1
        while len(a)>d:
            lead = a.pop()
            shift = len(a)-d
            for j in range(d):
                a[shift+j] -= lead*self.modulus[j]
        return tuple(a+[F(0)]*(d-len(a)))

    def add(self,a,b):
        return tuple(x+y for x,y in zip(a,b))

    def mul(self,a,b):
        raw = [F(0)]*(len(a)+len(b)-1)
        for i,x in enumerate(a):
            for j,y in enumerate(b):
                raw[i+j] += x*y
        # norm accepts Fraction internally only via exact string conversion.
        return self.norm([str(x) for x in raw])

    def inverse(self,a):
        if a == self.zero:
            raise ValueError('Division by zero')
        q = sp.Poly(sum(sp.Rational(str(c))*self.x**i for i,c in enumerate(a)),self.x,domain=sp.QQ)
        inv = sp.invert(q,self.p)
        return self.norm([str(c) for c in reversed(inv.all_coeffs())])

    def sign(self,a):
        if a == self.zero:
            return 0
        lo,hi = self.lo,self.hi
        for _ in range(4096):
            lower=upper=F(0)
            for c in reversed(a):
                products=[lower*lo,lower*hi,upper*lo,upper*hi]
                lower,upper=min(products)+c,max(products)+c
            if lower>0:
                return 1
            if upper<0:
                return -1
            mid=(lo+hi)/2
            if self.p.eval(mid)==0:
                value=sum(c*mid**i for i,c in enumerate(a))
                return (value>0)-(value<0)
            if self.p.count_roots(lo,mid)==1:
                hi=mid
            else:
                lo=mid
        raise ArithmeticError('Sign isolation exceeded limit; unknown, not accepted')


def verify(data):
    field=RealField(data['minpoly'],data['interval'])
    A=[field.norm(c) for c in data['A']]
    B=[field.norm(c) for c in data['B']]
    if not A or not B or A[0]!=field.one or B[0]!=field.one:
        raise ValueError('Constant coefficients must both be one')
    if A[-1]!=field.one or B[-1]!=field.one:
        raise ValueError('Leading coefficients must both be one')
    fractional=False
    for a in A+B:
        s=field.sign(a)
        below=field.sign(tuple(x-y for x,y in zip(field.one,a)))
        if s<0 or below<0:
            raise ValueError('Coefficient is outside [0,1]')
        fractional |= s>0 and below>0
    C=[field.zero]*(len(A)+len(B)-1)
    for i,a in enumerate(A):
        for j,b in enumerate(B):
            C[i+j]=field.add(C[i+j],field.mul(a,b))
    if any(c not in (field.zero,field.one) for c in C):
        raise ValueError('Product is not exactly a 0-1 polynomial')
    support=[i for i,c in enumerate(C) if c==field.one]
    sums=[tuple(sum(v[i] for v in P) for i in range(len(field.one))) for P in (A,B)]
    probabilities=[]
    for P,s in zip((A,B),sums):
        inv=field.inverse(s)
        probabilities.append([[str(v) for v in field.mul(a,inv)] for a in P])
    # Independently check the probability product identity inside the same field.
    assert field.mul(sums[0],sums[1])==field.norm([len(support)])
    return {'valid_factorization':True,'counterexample':bool(fractional),
            'sum_support':support,'uniform_sum_probability':str(F(1,len(support))),
            'probability_coefficients_in_theta':probabilities}


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('certificate',type=Path)
    args=ap.parse_args()
    print(json.dumps(verify(json.loads(args.certificate.read_text())),indent=2))
