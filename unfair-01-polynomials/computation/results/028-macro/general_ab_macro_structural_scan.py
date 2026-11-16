#!/usr/bin/env python3
"""Small exact scan for the conditional a<b<2a macro lemma.

It searches only the displayed finite box.  Its role is to find whether the
gcd hypothesis is automatic under the support/degree hypotheses, not to make
an unbounded claim.
"""
from collections import Counter, defaultdict
import sympy as sp
from itertools import combinations


def P(a):
    return Counter({(0,0):1,(0,a):1,(0,2*a):1,(a,-a):1,(a,0):1})


def mul(A,B):
    C=Counter()
    for (r,e),m in A.items():
        for (s,f),n in B.items(): C[r+s,e+f]+=m*n
    return C


def audit(a,b,d):
    Tcnt=mul(P(a),P(b)); T=set(Tcnt); TT=set(mul(Tcnt,Tcnt))
    if set(Tcnt.values())!={1}: return None
    im=defaultdict(list)
    for r,e in T: im[r*d+e].append((r,e))
    if not all(len(q)==1 for q in im.values()): return None
    aliases=[]
    for r,e in TT:
        if r*d+e in im and im[r*d+e][0]!=(r,e): aliases.append(((r,e),im[r*d+e][0]))
    if not {(x[0][0],x[1][0]) for x in aliases} <= {(a,b),(2*a,a+b)}: return None
    return aliases


u=sp.symbols('u')
# Formal (a,b)-independent row syzygy used in the conditional proof.
La,Lb,Ua,Ub,D,E,z=sp.symbols('La Lb Ua Ub D E z')
V=La*Lb; R=Ua*Lb+z*D; S=La*Ub-D; T=z*E; W=Ua*Ub-E
assert sp.expand(R*S*W-S**2*T-W**2*V-(La*E-Ua*D)*(Lb*W-z*Ub*S))==0
found=[]; total=0; gcd_fail=0; genuine=[]
for a in range(2,16):
    for b in range(a+1,2*a):
        for d in range(2,101):
            h=(b-a)*d
            aliases=audit(a,b,d)
            if aliases is None: continue
            total+=1
            La=1+u**a+u**(2*a); Lb=1+u**b+u**(2*b)
            H=lambda q: 1+q+q**2+q**(d-1)+q**d
            g1=sp.gcd(La,H(u**b)); g2=sp.gcd(Lb,H(u**a))
            if g1!=1 or g2!=1:
                gcd_fail+=1
                if len(found)<10: found.append((a,b,d,h,sp.factor(g1),sp.factor(g2),aliases))
                g=sp.gcd(La*Lb,sp.expand(H(u**a)*H(u**b)))
                irreducibles=[p for p,e in sp.factor_list(g)[1] for _ in range(e)]
                for mask in range(1,1<<len(irreducibles)):
                    q=sp.prod(p for j,p in enumerate(irreducibles) if mask>>j & 1)
                    qq=sp.Poly(q,u); coeff=qq.all_coeffs()
                    if not all(c>=0 for c in coeff): continue
                    quotient,rem=sp.div(sp.expand(H(u**a)*H(u**b)),q,u)
                    if rem==0 and all(c>=0 for c in sp.Poly(quotient,u).all_coeffs()):
                        genuine.append((a,b,d,sp.factor(q),aliases))
                        break
print('tested macro-admissible triples (no h threshold):',total)
print('gcd failures:',gcd_fail)
print('genuine nonnegative pure divisors:',genuine[:10])
print('general row syzygy verified symbolically')
for x in found: print(x)
