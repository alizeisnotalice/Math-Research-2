"""Nontrivial true Boolean selector spectrum, compressed exactly by setting address-core variables to zero."""
from collections import defaultdict, Counter
from fractions import Fraction
import json,math,platform

p=5
N=1<<p
h=13
d=h+p+N
C=1<<h
s=C+2*N*N
m=Fraction(s,d)
k=math.ceil(math.log2(float(m)*math.log2(float(m))))
M=1<<k
assert k<h
# N independent tail labels, p independent address directions, h core directions.
# The first core bit is the selector t; Walsh coefficients below have common denominator C.
S={}
for a in range(C):
    S[a]=(2-C//2) if a>>1==0 else 2
for i in range(N):
    tag=1<<(h+p+i)
    for addr in range(N):
        for t in range(2):
            alpha=tag|(addr<<h)|t
            S[alpha]=(C//(2*N))*(1-2*((t+(i&addr).bit_count())%2))
assert len(S)==s and sum(a*a for a in S.values())==C*C
core_merit=C-2*m*h
for j in range(h+1):
    assert (1<<j)-2*m*j<=core_merit
assert 2*N<=2*m
# For arbitrary T with j=dim(T intersection C), e outside label-projection rank:
# |T cap S|<=2^j+2Ne, dim T>=j+e, so these checks certify the global optimizer C.
new=defaultdict(int)
for alpha,coeff in S.items():
    beta=(alpha&((1<<k)-1))|((alpha>>h)<<k)
    new[beta]+=coeff
new={a:v for a,v in new.items() if v}
assert sum(x*x for x in new.values())==C*C
def rank(xs):
    piv={}
    for x in xs:
        while x.bit_length() in piv:
            x^=piv[x.bit_length()]
        if x:
            piv[x.bit_length()]=x
    return len(piv)
assert rank(S)==d and rank(new)==k+p+N
classes=0
for j in range(k+1):
    for e in range(N+1):
        assert (1<<j)-1+2*N*e<=2*m*(j+e)
        classes+=1
print(json.dumps(dict(python=platform.python_version(),address_bits=p,tail_labels=N,core_dim=h,
  proper_core_selected=True,core_merit=str(core_merit),s=s,d=d,compressed_kernel_dim=k,M=M,
  s_prime=len(new),d_prime=rank(new),analytic_dimension_class_checks=classes,
  original_actually_a_record=rank(new)**2*s<=d*d*len(new),
  check="exact normalized Walsh via common denominator, no enumeration of the 50-dimensional lattice"),indent=2))
