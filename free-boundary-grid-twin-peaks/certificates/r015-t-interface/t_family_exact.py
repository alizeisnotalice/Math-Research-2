#!/usr/bin/env python3
"""Finite exact T-family masses and root moments, using subset integers.

Only n=3,5,7.  The all-n exclusion is analytic, not a trend extrapolation.
"""
from fractions import Fraction
from math import comb, factorial
from pathlib import Path
import json
import sys


def rectangle(w, h):
    N = w*h
    adj = [0]*N
    for v in range(N):
        x,y = v%w,v//w
        for dx,dy in ((1,0),(-1,0),(0,1),(0,-1)):
            X,Y=x+dx,y+dy
            if 0<=X<w and 0<=Y<h: adj[v] |= 1 << (Y*w+X)
    a,mx,my,m2 = ([0]*(1<<N) for _ in range(4))
    for v in range(N):
        bit=1<<v;x,y=v%w,v//w
        a[bit]=1;mx[bit]=x;my[bit]=y;m2[bit]=x*x+y*y
    for S in range(1,1<<N):
        if S&(S-1)==0: continue
        bits=S
        while bits:
            b=bits&-bits;bits-=b;v=b.bit_length()-1;T=S^b
            if adj[v]&T:
                a[S]+=a[T];mx[S]+=mx[T];my[S]+=my[T];m2[S]+=m2[T]
    return a[-1],mx[-1],my[-1],m2[-1]


def main():
    # These total A2 integers are independently certified in rounds 1 and 5.
    A2={3:143808,5:7063261481362958280000,
        7:137077928808481869759712347602605837601908319657172480}
    rows=[]
    for n in (3,5,7):
        m=(n-1)//2;a=2*m*m;g=4*m+1;k=g-2
        count,mx,my,m2=rectangle(m,2*m)
        assert Fraction(mx,count)==Fraction(m-1,2)
        assert Fraction(my,count)==Fraction(2*m-1,2)
        var=Fraction(m2,count)-Fraction(mx,count)**2-Fraction(my,count)**2
        suffix=Fraction(6*m+1,4)*factorial(g-1)
        assert suffix.denominator==1
        family=count**2*comb(2*a,a)*suffix.numerator
        mass=Fraction(family,A2[n]);assert 0<mass<=1
        q=(2*var+(m+1)**2)/n**2
        Ek=factorial(2*k)//(2**k*factorial(k))
        exact_new_bound=Fraction((2*a)*(2*a-1),2*a*a)*suffix/(2*Ek)
        assert mass<=exact_new_bound
        rows.append({'n':n,'rectangle':[m,2*m],'onepeak_rectangle':str(count),
                     'branch_root_variance':str(var),'T_family_count':str(family),
                     'A2_total':str(A2[n]),'T_family_mass':str(mass),
                     'T_family_conditional_Q':str(q),
                     'new_comparison_ratio_exact':str(exact_new_bound),
                     'mass_decimal_display':float(mass),'Q_decimal_display':float(q)})
    result={'status':'passed','arithmetic':'Python integers/Fraction; float display only','rows':rows}
    output=Path(sys.argv[1]) if len(sys.argv)>1 else Path('t-family-exact.json')
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
