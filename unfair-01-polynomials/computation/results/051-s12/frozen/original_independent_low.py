"""Independent low-row replay: direct rational-surd convolutions.

Uses no original low-screen functions or fixed integer windows. Enumerates
full O_K coefficients, derives free integer row choices from unit intervals,
then computes C2/C3 directly. All rejected lattice pairs are recomputed.
"""
import json,math,hashlib,time
from fractions import Fraction as F
from pathlib import Path
from collections import Counter
T=Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials')
def add(x,y):return(x[0]+y[0],x[1]+y[1])
def neg(x):return(-x[0],-x[1])
def mul(x,y,D):return(x[0]*y[0]+D*x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def scale(n,x):return(n*x[0],n*x[1])
def bar(x):return(x[0],-x[1])
def omega(m,n,D):return(F(2*m+n,2),F(n,2)) if D%4==1 else(F(m),F(n))
def sg(x,D):
    p,q=x
    if q==0:return(p>0)-(p<0)
    if p>=0 and q>0:return 1
    if p<=0 and q<0:return -1
    z=p*p-D*q*q;assert z!=0
    return ((z>0)-(z<0))*(1 if p>0 else -1)
def unit(x,D):return sg(x,D)>=0 and sg(add((1,0),neg(x)),D)>=0
def fl(x,D):
    lo,hi=-1,1
    while sg(add(x,(-lo,0)),D)<0:hi=lo;lo*=2
    while sg(add(x,(-hi,0)),D)>=0:lo=hi;hi*=2
    while hi-lo>1:
        mid=(lo+hi)//2
        if sg(add(x,(-mid,0)),D)>=0:lo=mid
        else:hi=mid
    return lo
def shifts(x,D):
    n=-fl(x,D)
    return [j for j in (n-1,n,n+1) if unit(add(x,(j,0)),D)]
def bound(x,D,b):return abs(x[0])+abs(x[1])*0<b and (b-abs(x[0]))>0 and x[1]*x[1]*D<(b-abs(x[0]))**2
def lattice(D,B,positive):
    M=math.ceil(2*B if D%4==1 else B)
    return [(m,n,omega(m,n,D)) for n in range(1 if positive else -M,M+1) if n for m in range(-M,M+1) if bound(omega(m,n,D),D,B)]
def coefficient(A,B,k,D):
    z=(F(0),F(0))
    for i in range(k+1):
        if i<len(A) and k-i<len(B):z=add(z,mul(A[i],B[k-i],D))
    return z
def main():
    start=time.monotonic();counts=Counter();rows=[]
    for D in range(2,289):
        if any(D%(i*i)==0 for i in range(2,math.isqrt(D)+1)):continue
        aa=lattice(D,F(17,2),True)
        if not aa:continue
        bb=lattice(D,F(117,4),False);counts['fields']+=1
        for am,an,a in aa:
            g1=-fl(a,D);A1=add(a,(g1,0));q1=1-2*a[0]-g1;assert q1.denominator==1;q1=int(q1);B1=add(bar(a),(q1,0))
            for bm,bn,b in bb:
                counts['raw']+=1
                if not bound(add(mul(a,a,D),scale(-2,b)),D,F(55,4)):counts['s2']+=1;continue
                A2base=add(b,scale(g1,a));B2base=add(bar(b),scale(q1,bar(a)))
                C2base=add(add(A2base,B2base),mul(A1,B1,D));assert C2base[1]==0 and C2base[0].denominator==1
                choices=[]
                for g2 in shifts(A2base,D):
                    q2=1-int(C2base[0])-g2;A2=add(A2base,(g2,0));B2=add(B2base,(q2,0))
                    if unit(B2,D):choices.append((g2,q2,A2,B2))
                assert len(choices)<=1
                if not choices:counts['A2']+=1;continue
                g2,q2,A2,B2=choices[0];A3base=add(scale(g1,b),scale(g2,a));B3base=add(scale(q1,bar(b)),scale(q2,bar(a)))
                cross=add(mul(A1,B2,D),mul(A2,B1,D));C3base=add(add(A3base,B3base),cross)
                assert C3base[1]==0 and C3base[0].denominator==1
                choices=[]
                for rho in shifts(A3base,D):
                    q0=1-int(C3base[0])-rho;A3=add(A3base,(rho,0));B3=add(B3base,(q0,0))
                    if unit(B3,D):choices.append((rho,q0))
                assert len(choices)<=1
                if not choices:counts['A3']+=1;continue
                rho,q0=choices[0];h=math.gcd(an,abs(bn));alpha,beta=an//h,bn//h
                if (rho-q0)%alpha:counts['high_congruence']+=1;continue
                counts['survivor']+=1;rows.append((D,am,an,bm,bn,g1,q1,g2,q2,rho,q0,alpha,beta))
    source=T/'agents/algebraic_global/round051-s12/ab-lowprefix-floor-exact/full/summary.json';x=json.loads(source.read_text());keys=('D','a_m','a_n','b_m','b_n','g1','q1','g2','q2','rho','q0','alpha','beta');expected=[tuple(z[k] for k in keys) for z in x['rows']]
    assert len(rows)==len(set(rows))==7245 and sorted(rows)==sorted(expected) and counts['raw']==396738
    out={'status':'passed_independent_all_lattice_pairs_and_7245_low_prefixes','counts':dict(counts),'expected_input_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'seconds':time.monotonic()-start}
    Path(__file__).with_name('independent-s12-low-audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
if __name__=='__main__':main()
