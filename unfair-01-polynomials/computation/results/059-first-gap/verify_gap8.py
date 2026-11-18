"""Exact elimination and independent finite-field irreducibility tests for C8."""
import argparse
import json
from math import gcd, isqrt
from pathlib import Path
import sympy as sp

t,z=sp.symbols('t z')
COEFFS={
    'F17':[19,-168,768,-2378,5488,-9738,13376,-14768,15044,-17688,23856,-29176,27664,-18928,8832,-2632,472,-48],
    'F10':[3,-24,92,-218,352,-404,328,-182,64,-14,2],
    'F5':[5,-20,40,-44,26,-6],
    'F6':[3,-6,8,-9,9,-5,2],
}
PRIMES={'F17':37,'F10':13,'F5':7,'F6':5}


def trim(a):
    while a and not a[-1]:
        a.pop()
    return a


def remainder(u,v,p):
    a=list(u)
    inv=pow(v[-1],-1,p)
    while a and len(a)>=len(v):
        k=len(a)-len(v)
        factor=a[-1]*inv%p
        for j,b in enumerate(v):
            a[k+j]=(a[k+j]-factor*b)%p
        trim(a)
    return a


def multiply(u,v,modulus,p):
    c=[0]*(len(u)+len(v)-1)
    for i,a in enumerate(u):
        for j,b in enumerate(v):
            c[i+j]=(c[i+j]+a*b)%p
    return remainder(trim(c),modulus,p)


def power(u,e,modulus,p):
    out=[1]
    while e:
        if e%2:
            out=multiply(out,u,modulus,p)
        u=multiply(u,u,modulus,p)
        e//=2
    return out


def gcd_poly(a,b,p):
    while b:
        a,b=b,remainder(a,b,p)
    return [(v*pow(a[-1],-1,p))%p for v in a]


def irreducible_certificate(coefficients,prime):
    assert all(prime%d for d in range(2,isqrt(prime)+1))
    assert gcd(*coefficients)==1
    f=[c%prime for c in coefficients[::-1]]
    assert f[-1]!=0
    n=len(f)-1
    divisors=[k for k in range(2,n+1) if n%k==0 and all(k%d for d in range(2,isqrt(k)+1))]
    needed={n//q for q in divisors}
    checks={}
    value=[0,1]
    for i in range(1,n+1):
        value=power(value,prime,f,prime)
        if i in needed:
            delta=value+[0]*max(0,2-len(value))
            delta[1]=(delta[1]-1)%prime
            g=gcd_poly(f,trim(delta),prime)
            assert g==[1]
            checks[i]=g
    assert value==[0,1]
    return {'prime':prime,'degree':n,'degree_prime_divisors':divisors,
            'gcd_x_p_power_minus_x':checks,'frobenius_final':value}


def main():
    pa=argparse.ArgumentParser()
    pa.add_argument('--out',type=Path,required=True)
    a=pa.parse_args()
    assert __debug__
    F={key:sum(c*t**(len(row)-1-i) for i,c in enumerate(row)) for key,row in COEFFS.items()}
    m=1-t+t*t
    n=1-t*m+(2*t-1)*z
    p=1-t*n-z*(m-z)
    q=1-(1-t)*p-z*n
    eq6=t*q+p*(m-z)-1
    eq7=z*q+p*n-1
    assert sp.expand(sp.resultant(eq6,eq7,z)+t**4*F['F17'])==0
    J=t**5-3*t**4+2*t**3-t**2+t-1
    for side in (0,1):
        u=m if side==0 else 0
        w=m-u
        n=1-t*w-(1-t)*u
        p=1-t*(n-z)-(1-t)*z-u*w
        r=1-p*w-z*(n-z)
        eq5=(1-t)*p+u*(n-z)+w*z-1
        eq7=((1-t) if side==0 else t)*r+p*(n-z)-1
        result=t*F['F10'] if side==0 else t**2*(t**2+2)*(t**2-2*t+2)*J
        assert sp.expand(sp.resultant(eq5,eq7,z)-result)==0
    assert sp.expand(J-t**3*(t-1)*(t-2)+m)==0
    pure={(0,0):-t*((1-t)**2+1)*((1-t)**4+1),
          (0,1):-t**2*F['F5'],
          (1,0):(t-1)*(t**2+1)*(t**4+1),
          (1,1):-t*F['F6']}
    for (i,j),result in pure.items():
        A=[1,t];B=[1,1-t]
        for k,side in enumerate((i,j,0,j,i),2):
            free=sp.expand(1-sum(A[h]*B[k-h] for h in range(1,k)))
            A.append(free if side==0 else 0)
            B.append(free if side==1 else 0)
        eq7=sp.expand(sum(A[h]*B[7-h] for h in range(1,7))-1)
        assert sp.expand(eq7-result)==0
    output={'status':'verified_all_gap8_eliminants','sympy':sp.__version__,
            'irreducibility':{key:irreducible_certificate(row,PRIMES[key]) for key,row in COEFFS.items()},
            'primitive_nonmonic_coefficients_descending':COEFFS,
            'mixed_resultants_checked':3,'single_side_patterns_checked':4}
    with a.out.open('x') as f:
        f.write(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output,sort_keys=True))


if __name__=='__main__':
    main()
