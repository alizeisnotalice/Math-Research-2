#!/usr/bin/env python3
"""Sparse exact verification for the p=2,q=17 coprime-slice example."""
from __future__ import annotations
from functools import reduce
from math import gcd
import hashlib, json, os, sys
from pathlib import Path

assert os.path.ismount('/Volumes/data')

def mul(a, b):
    z = {}
    for i,u in a.items():
        for j,v in b.items(): z[i+j] = z.get(i+j,0)+u*v
    return {i:v for i,v in z.items() if v}

def subst(a, d): return {d*i:v for i,v in a.items()}
def plus(a, b):
    z=dict(a)
    for i,v in b.items(): z[i]=z.get(i,0)+v
    return {i:v for i,v in z.items() if v}
def degree(a): return max(a)
def binary(a): return set(a.values()) == {1}
def section(a, q):
    z={}
    for i,v in a.items():
        if i%q==0: z[i//q]=z.get(i//q,0)+v
    return {i:v for i,v in z.items() if v}

# Polynomials over F_2 are nonnegative integers whose bits are coefficients.
def degree2(a): return a.bit_length()-1
def rem2(a, f):
    while a and degree2(a)>=degree2(f): a ^= f << (degree2(a)-degree2(f))
    return a
def mul2(a,b):
    z=0
    while b:
        if b&1:z^=a
        a<<=1;b>>=1
    return z
def gcd2(a,b):
    while b: a,b=b,rem2(a,b)
    return a
def subst2(a,k):
    z=0; i=0
    while a:
        if a&1:z ^= 1 << (i*k)
        a >>= 1; i += 1
    return z
def pow2_raw(a,n):
    z=1
    while n:
        if n&1:z=mul2(z,a)
        a=mul2(a,a);n>>=1
    return z
def bits_from_sparse(a):
    z=0
    for i,c in a.items():
        if c&1:z ^= 1 << i
    return z
def mulmod2(a,b,f): return rem2(mul2(a,b),f)
def powmod2(a,n,f):
    z=1
    while n:
        if n&1:z=mulmod2(z,a,f)
        a=mulmod2(a,a,f);n>>=1
    return z
def p_of_power(k, f):
    z=powmod2(2,k,f)  # x^k in F_2[x]/(f)
    return 1 ^ z ^ mulmod2(mulmod2(z,z,f),z,f)

def main():
    p,q,N=2,17,5001
    P={0:1,1:1,3:1}
    S0=mul(P,subst(P,4))
    F=mul(subst(P,17),{0:1,N:1})
    G=subst(P,68)
    H=mul(F,G)
    RF,RG=section(F,q),section(G,q)
    U=mul(F,subst(RG,p)); V=mul(G,subst(RF,p))
    A=mul(subst(F,p),subst(G,q));B=mul(subst(F,q),subst(G,p));C=mul(subst(H,p),subst(H,q))
    assert H==mul(subst(S0,17),{0:1,N:1})
    assert RF==P and RG==subst(P,4) and S0==mul(RF,RG)
    assert all(binary(z) for z in [H,U,V,A,B,C])
    assert C==mul(A,B) and mul(U,V)==mul(H,subst(S0,p)) and binary(mul(U,V))
    assert degree(F)==5052 and degree(G)==204 and min(degree(F),degree(G))>=2*q
    assert reduce(gcd, set(F)|set(G))==1
    baseC=mul(subst(S0,34),subst(S0,289)); assert binary(baseC) and degree(baseC)==4845
    assert min(2*N,17*N-2*N)>degree(baseC)
    # Mod 2: P=x^3+x+1 is irreducible because it has neither F2 root.
    f=0b1011
    assert rem2(f,0b10)==1 and rem2(f,0b11)==1
    # Directly construct S(x^2)=P(x^2)P(x^8)=P(x)^10, then compute the
    # Euclidean gcd with the sparse H reduction.  This avoids inferring the
    # last factor of H from P(alpha^N): that factor is 1+alpha^N instead.
    sxp_bits=mul2(subst2(f,2),subst2(f,8))
    assert sxp_bits==pow2_raw(f,10)
    hbits=bits_from_sparse(H)
    assert gcd2(hbits,sxp_bits)==1
    P_at_17,P_at_68=p_of_power(17,f),p_of_power(68,f)
    alpha_N=powmod2(2,N,f)
    hmod=rem2(hbits,f)
    factor_hmod=mulmod2(mulmod2(P_at_17,P_at_68,f),1^alpha_N,f)
    assert P_at_17!=0 and P_at_68!=0 and (1^alpha_N)!=0
    assert hmod==factor_hmod and hmod!=0
    out=dict(status='passed_sparse_coprime_slice_descent_example',
             parameters=dict(p=p,q=q,N=N),
             degrees=dict(F=degree(F),G=degree(G),H=degree(H),A=degree(A),B=degree(B),C=degree(C),U=degree(U),V=degree(V),base_C=degree(baseC)),
             support_sizes={name:len(poly) for name,poly in dict(F=F,G=G,H=H,A=A,B=B,C=C,U=U,V=V,S=S0).items()},
             binary_checked=['H','A','B','C','U','V','UV'],
             sections=dict(RF=sorted(RF),RG=sorted(RG),S=sorted(S0)),
             joint_support_gcd=reduce(gcd,set(F)|set(G)),
             mod2=dict(P_irreducible_degree_3=True,order_of_root=7,exponents_mod_7=dict(x17=3,x68=5,N=3),P_at_alpha17=P_at_17,P_at_alpha68=P_at_68,alpha_to_N=alpha_N,one_plus_alpha_to_N=1^alpha_N,H_mod_P=hmod,direct_f2_gcd_degree=degree2(gcd2(hbits,sxp_bits)),Sxp_equals_P_to_10=True,gcd_H_Sxp_over_Q=True),
             scope='Exact verification of one full crossed, double-long example satisfying the added coprimality hypothesis; not a proof that the hypothesis is automatic.',
             interpreter=sys.version,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps(out,sort_keys=True))

if __name__=='__main__': main()
