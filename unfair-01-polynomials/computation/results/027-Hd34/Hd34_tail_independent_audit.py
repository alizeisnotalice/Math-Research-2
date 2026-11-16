#!/usr/bin/env python3
"""Independent exact audit for the H_d(x^3)H_d(x^4), d>=17 macro tail.

No external source code or JSON is imported.  The finite set operations certify
the support/carry accounting only; the accompanying review supplies the
general nonnegative-factor proof.
"""
from collections import Counter, defaultdict
import sympy as sp


def p(a):
    return Counter({(0,0):1, (0,a):1, (0,2*a):1, (a,-a):1, (a,0):1})


def multiply(a,b):
    c=Counter()
    for (r,e),m in a.items():
        for (s,f),n in b.items(): c[r+s,e+f]+=m*n
    return c


Tcounter=multiply(p(3),p(4)); T=set(Tcounter)
assert sum(Tcounter.values())==25 and set(Tcounter.values())=={1}
assert {r:sorted(e for rr,e in T if rr==r) for r in {0,3,4,7}} == {
    0:[0,3,4,6,7,8,10,11,14], 3:[-3,0,1,4,5,8],
    4:[-4,-1,0,2,3,6], 7:[-7,-4,-3,0]}
TT=set(multiply(Tcounter,Tcounter))


def image(d):
    out=defaultdict(list)
    for r,e in T: out[r*d+e].append((r,e))
    return out


def aliases(d):
    im=image(d); assert all(len(v)==1 for v in im.values()), (d,im)
    return sorted(((r,e),im[r*d+e][0]) for r,e in TT
                  if r*d+e in im and (r,e)!=im[r*d+e][0])


def all_parameters(left,right):
    vals=set()
    for r,e in left:
        for s,f in right:
            if r!=s and (f-e)%(r-s)==0 and (f-e)//(r-s)>0:
                vals.add((f-e)//(r-s))
    return vals


assert {d for d in all_parameters(T,T) if d>=11}=={12}
par=all_parameters(TT,T)
assert max(par)==26 and not {d for d in par if d>=27}
need={14:{-4},15:{-7,-3},17:{-7,-4},18:set(),19:{-7,-3},20:{-7,-4},21:set(),22:set(),23:{-7},24:set(),25:set(),26:set()}
for d,want in need.items():
    got=aliases(d)
    assert {(a[0],c[0]) for a,c in got} <= {(3,4),(6,7)},(d,got)
    supplied={c[1] for a,c in got if a[0]==6 and c[0]==7}
    assert supplied==want,(d,supplied,want)
    assert supplied != {-7,-4,-3,0}
    # For nonzero S, deg(u^d U4 S) >= d-4; its degree exceeds <=8 from L4 W.
    assert d-4 > 8

# The remaining small macro patterns relevant to the strengthened argument.
for d, want7, want4 in [
    (11, {-7,-3}, {-4,-1,0,2,3}),
    (13, {-7,-4,-3,0}, {-4,-1,0,2,3,6}),
    (16, {-7,-4,-3,0}, {-4,-1,0,2,3,6}),
]:
    got=aliases(d)
    r7={c[1] for a,c in got if a[0]==6 and c[0]==7}
    r4={c[1] for a,c in got if a[0]==3 and c[0]==4}
    assert r7==want7 and r4==want4,(d,r7,r4)
    assert {(a[0],c[0]) for a,c in got} <= {(3,4),(6,7)}

# Exact formal row identity.
L3,L4,U3,U4,D,E,z=sp.symbols('L3 L4 U3 U4 D E z')
V=L3*L4; R=U3*L4+z*D; S=L3*U4-D; Trow=z*E; W=U3*U4-E
assert sp.expand(R*S*W-S**2*Trow-W**2*V-(L3*E-U3*D)*(L4*W-z*U4*S))==0

# Exact gcd checks: d-period 36 covers all root-of-unity factors of L3*L4.
u=sp.symbols('u'); ell3=1+u**3+u**6; ell4=1+u**4+u**8
for residue in range(36):
    d=residue or 36
    h=lambda t:1+t+t**2+t**(d-1)+t**d
    assert sp.gcd(ell3*ell4,sp.expand(h(u**3)*h(u**4)))==1
    assert sp.gcd(ell4,h(u**4))==1

print('all actual-to-target carry parameters <=26; target collision >=11 only d=12')
for d in need: print(d, sorted(need[d]))
print('d=11 r7/r4 carries:', sorted({c[1] for a,c in aliases(11) if a[0]==6}),
      sorted({c[1] for a,c in aliases(11) if a[0]==3}))
print('d=13,16 have full r7 and full r4 carry; this is the W=0 exceptional pattern')
print('row syzygy and 36-period gcds passed exactly')
