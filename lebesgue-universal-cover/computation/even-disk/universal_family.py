from fractions import Fraction as F
from flint import arb,ctx
import json,sys,platform
ctx.prec=192;M=10**6;P=[634599,971435,997893,999779,1096241,1007040,1000425]
EXPECTED={(634599,1096241,461642),(634599,1096241,1730840),(971435,1007040,958889),(971435,1096241,1144392),(997893,1000425,997026),(997893,1007040,1010439),(999779,1000000,999779),(999779,1000425,1000646),(1000000,1000000,1000000)}
def fingerprint(n):
 m=n//2;w=[M]*n;w[0]=w[1]=P[0]
 for k in range(3):w[2+k]=w[n-1-k]=P[1+k];w[m-k]=w[m+1+k]=P[4+k]
 W=sum(w);prefix=[0]
 for z in w:prefix.append(prefix[-1]+z)
 b=[2*prefix[i]+w[i] for i in range(n)];B=b+[x+2*W for x in b]
 cases=[('low',q,q) for q in range(10)]+[('mid',q,m+q) for q in range(-10,10)]+[('high',q,n+q) for q in range(-10,0)]+[('leftbulk',0,10),('rightbulk',0,m+10)]
 return {(name,q,k):(w[i],w[(i+m+k)%n],B[i+m+k]-b[i]-W) for name,q,i in cases for k in range(-2,3)}
assert fingerprint(2000)==fingerprint(2002)
def check(n):
 assert n%2==0 and n>=2000;m=n//2;w=[M]*n;w[0]=w[1]=P[0]
 for k in range(3):w[2+k]=w[n-1-k]=P[1+k];w[m-k]=w[m+1+k]=P[4+k]
 W=sum(w);assert W==M*n-585176;assert sum(abs(x-M) for x in w)==M
 prefix=[0]
 for z in w:prefix.append(prefix[-1]+z)
 b=[2*prefix[i]+w[i] for i in range(n)];B=b+[x+2*W for x in b];types=set();cases=set(range(10))|set(range(m-10,m+10))|set(range(n-10,n))|{10,m+10}
 for i in cases:
  for k in range(-2,3):
   j=i+m+k;u=abs(B[j]-b[i]-W);p=tuple(sorted((w[i],w[j%n])))+(u,)
   if u<2*M:types.add(p)
 assert types==EXPECTED,(n,types^EXPECTED)
 # Eight strict type bounds and one equality, purely integer.
 gaps=[3*M*M-(2*(a*a+b*b)-d*d) for a,b,d in sorted(types) if (a,b,d)!=(M,M,M)];assert min(gaps)==441951159
 # t=pi/(n-.585176) < 1/500; coefficient remainder bounded by 51t².
 assert F(441951159,M*M)>F(51,500**2)
 assert 6*M-(2*M+(max(P)-min(P))+585176)>2*M
 t=arb.pi()/(n-arb(73147)/125000);assert t<arb(1)/500
 D=2/(t.cos())*(t/2).cos();e=1/D*2
 return {'n':n,'W':W,'representative_indices':len(cases),'type_count':len(types),'minimum_integer_quadratic_gap':min(gaps),'comparison':str(D),'incircle_diameter_ratio':str(e),'tail_spacing_proof_limit':'t<1/500'}
rows=[check(n) for n in (2000,2002,4302,6084,7790,10000)];out={'precision_bits':ctx.prec,'python':platform.python_version(),'scope':'fixed integer-weight family; proof of uniform type list and inequalities is in README','symbolic_representative_fingerprints_equal_for_2000_2002':True,'alpha':'73147/125000','exact_third_order_upper_coefficient_over_pi_squared':'219441/500000','near_pair_types':[list(x) for x in sorted(EXPECTED)],'checks':rows};open(sys.argv[1],'w').write(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
