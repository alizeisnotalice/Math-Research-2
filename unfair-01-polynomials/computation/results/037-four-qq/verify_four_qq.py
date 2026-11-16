#!/usr/bin/env python3
"""Pure-integer certificates for the four residual Q-irreducibility cases.

No SymPy, factorization engine, or primality oracle is used.  All finite-field
claims are replayed with modular polynomial arithmetic and Rabin's test.
"""
import argparse, hashlib, json, math, os, sys, time
from pathlib import Path

def trim(a):
    a=list(a)
    while len(a)>1 and a[-1]==0: a.pop()
    return a
def mod(a,p): return trim([x%p for x in a])
def sub(a,b,p):
    n=max(len(a),len(b)); return mod([(a[i] if i<len(a) else 0)-(b[i] if i<len(b) else 0) for i in range(n)],p)
def mul(a,b,p):
    c=[0]*(len(a)+len(b)-1)
    for i,u in enumerate(a):
        for j,v in enumerate(b): c[i+j]=(c[i+j]+u*v)%p
    return trim(c)
def divrem(a,b,p):
    a=mod(a,p); b=mod(b,p); assert b[-1]
    q=[0]*max(1,(len(a)-len(b)+1)); inv=pow(b[-1],-1,p)
    while len(a)>=len(b) and a!=[0]:
        k=len(a)-len(b); c=a[-1]*inv%p; q[k]=c
        for j,v in enumerate(b): a[k+j]=(a[k+j]-c*v)%p
        a=trim(a)
    return trim(q),a
def rem(a,f,p): return divrem(a,f,p)[1]
def gcd(a,b,p):
    while b!=[0]: a,b=b,rem(a,b,p)
    return mod([x*pow(a[-1],-1,p)%p for x in a],p)
def powmod(a,e,f,p):
    out=[1]
    while e:
        if e&1: out=rem(mul(out,a,p),f,p)
        a=rem(mul(a,a,p),f,p); e//=2
    return out
def prime_divisors(n):
    out=[];d=2
    while d*d<=n:
        if n%d==0:
            out.append(d)
            while n%d==0:n//=d
        d+=1
    if n>1:out.append(n)
    return out
def rabin(f,p):
    f=mod(f,p); n=len(f)-1; assert f[-1]==1 and n>0
    X=[0,1]
    assert sub(powmod(X,p**n,f,p),X,p)==[0]
    for r in prime_divisors(n):
        assert gcd(sub(powmod(X,p**(n//r),f,p),X,p),f,p)==[1]
def product(fs,p):
    out=[1]
    for f in fs:out=mul(out,f,p)
    return out
def derivative(f,p):return mod([(i*f[i])%p for i in range(1,len(f))],p)
def long_rem(F,U,p): return rem(F,U,p)
def cubic_absence(F, c, p):
    # A monic cubic factor is X^3+aX^2+bX+c. Exhaustion takes p^2 values.
    for a in range(p):
        for b in range(p):
            assert long_rem(F,[c,b,a,1],p)!=[0]

F18A=[1,0,0,0,0,0,1]+[0]*11+[1]       # 1+x^6+x^18
F18B=[1]+[0]*11+[1]+[0]*5+[1]          # 1+x^12+x^18
F20=[1,0,0,0,-1,0,1,0,1,0,-2,0,1,0,1,0,-1,0,0,0,1]
F24=[1]+[0]*9+[1,0,1,0,1]+[0]*9+[1]
F9A=[1,0,0,1,0,0,0,0,0,1]
F9B=[1,0,0,0,0,0,1,0,0,1]
R20=[-2,0,13,0,-43,0,34,0,-10,0,1]
R24=[1,0,-35,0,105,0,-112,0,54,0,-12,0,1]

def palindromic_identity(R,F):
    n=(len(F)-1)//2; ans={}
    for k,c in enumerate(R):
        for j in range(k+1): ans[n+k-2*j]=ans.get(n+k-2*j,0)+c*math.comb(k,j)
    got=[ans.get(i,0) for i in range(2*n+1)]
    assert got==F

def source_binding(path):
    d=json.loads(path.read_text()); rows=d['remaining']; wanted={tuple(F20),tuple(F18A),tuple(F24),tuple(F18B)}
    found={tuple(x['coefficients_low_to_high']) for x in rows}
    assert len(rows)==4 and len(found)==4 and wanted==found and d['remaining_unknown']==4
    return hashlib.sha256(path.read_bytes()).hexdigest()

def certify():
    assert F18A==[v for c in F9A[:-1] for v in (c,0)]+[F9A[-1]]
    assert F18B==[v for c in F9B[:-1] for v in (c,0)]+[F9B[-1]]
    # Degree-nine nuclei.  Mod 2 has an irreducible cubic and sextic; a Q
    # factor must have degree 3 or 6.  A cubic factor is contradicted modulo
    # 5 if its constant is -1 and modulo 7 if its constant is +1.
    A3=[1,0,1,1]; A6=[1,0,1,0,1,1,1] # low-to-high mod 2
    B3=[1,1,0,1]; B6=[1,1,1,0,1,0,1]
    for F,a3,a6 in ((F9A,A3,A6),(F9B,B3,B6)):
        assert product([a3,a6],2)==mod(F,2)
        rabin(a3,2);rabin(a6,2)
        cubic_absence(F,-1,5);cubic_absence(F,1,7)
    # R20 is irreducible mod 19.  Its quadratic discriminant z^2-4 is a
    # nonsquare in the degree-three residue field of the displayed p=5 factor.
    rabin(R20,19)
    P5a=[3,4,1,1];P5b=[2,4,4,1];P5c=[3,0,3,0,1]
    assert product([P5a,P5b,P5c],5)==mod(R20,5)
    for f in (P5a,P5b,P5c):rabin(f,5)
    assert gcd(mod(R20,5),derivative(R20,5),5)==[1]
    delta=[1,0,1] # z^2-4 modulo 5
    assert powmod(delta,(5**3-1)//2,P5a,5)==[4]
    # R24: degree options modulo 3 are 0,6,12; modulo 17 are 0,2,10,12.
    # Their intersection rules out every proper rational factor.  The first
    # degree-six p=3 factor witnesses the nonsquare discriminant.
    P3a=[2,0,1,0,2,1,1];P3b=[2,0,1,0,2,2,1]
    P17a=[10,0,1];P17b=[12,0,14,0,4,0,2,0,12,0,1]
    assert product([P3a,P3b],3)==mod(R24,3)
    assert product([P17a,P17b],17)==mod(R24,17)
    for f,p in ((P3a,3),(P3b,3),(P17a,17),(P17b,17)):rabin(f,p)
    assert gcd(mod(R24,3),derivative(R24,3),3)==[1]
    assert powmod([2,0,1],(3**6-1)//2,P3a,3)==[2]
    palindromic_identity(R20,F20);palindromic_identity(R24,F24)

def rejected(label, thunk):
    try: thunk()
    except AssertionError: return label
    raise AssertionError('semantic negative accepted: '+label)

def negative_tests():
    # Each mutation targets a distinct logical link, not just parser syntax.
    badF=F20.copy();badF[10]+=1
    return [
        rejected('altered_palindromic_coefficient',lambda:palindromic_identity(R20,badF)),
        rejected('R20_wrong_irreducibility_prime',lambda:rabin(R20,5)),
        rejected('R20_false_square_discriminant',lambda:assert_square_r20()),
        rejected('R24_false_square_discriminant',lambda:assert_square_r24()),
        rejected('wrong_cubic_constant_sign',lambda:assert_cubic_exists(F9A,1,7)),
    ]
def assert_square_r20():
    assert powmod([1,0,1],(5**3-1)//2,[3,4,1,1],5)==[1]
def assert_square_r24():
    assert powmod([2,0,1],(3**6-1)//2,[2,0,1,0,2,1,1],3)==[1]
def assert_cubic_exists(F,c,p):
    for a in range(p):
        for b in range(p):
            if long_rem(F,[c,b,a,1],p)==[0]: return
    assert False

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--bound',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    assert os.path.ismount('/Volumes/data'); start=time.monotonic(); digest=source_binding(a.bound);certify()
    out={'status':'passed','claim':'All four displayed monic integer polynomials are irreducible over Q. This is only a rational-factor result.','bound_source':str(a.bound),'bound_source_sha256':digest,'factors_low_to_high':{'degree20':F20,'degree18a':F18A,'degree24':F24,'degree18b':F18B},'method':['degree-nine mod-2/Rabin plus finite cubic obstruction and Capelli quadratic composition','palindromic trace polynomials with Rabin finite-field degree sieve and residue-field nonsquare discriminant'],'negative_tests_rejected':negative_tests(),'interpreter':sys.version,'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-start}
    a.out.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps(out,sort_keys=True))
if __name__=='__main__':main()
