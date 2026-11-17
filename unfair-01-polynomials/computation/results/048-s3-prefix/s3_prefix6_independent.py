"""Complete prefix<=6 witness enumeration via actual A/B rows, integer signs."""
import hashlib,json,math,time
from pathlib import Path

def add(*xs):return sum(x[0] for x in xs),sum(x[1] for x in xs)
def mul(x,y,T,N):return x[0]*y[0]-N*x[1]*y[1],x[0]*y[1]+x[1]*y[0]+T*x[1]*y[1]
def scale(x,n):return x[0]*n,x[1]*n
def sign(x,T,D):
    L=2*x[0]+T*x[1];q=x[1]
    if q==0:return (L>0)-(L<0)
    if L>=0 and q>0:return 1
    if L<=0 and q<0:return -1
    if L>0:return (L*L>q*q*D)-(L*L<q*q*D)
    return (q*q*D>L*L)-(q*q*D<L*L)
def floor(x,T,D):
    if x[1]==0:return x[0]
    lo,hi=-10000,10000
    assert sign(add(x,(-lo,0)),T,D)>0 and sign(add(x,(-hi,0)),T,D)<0
    while hi-lo>1:
        mid=(hi+lo)//2
        if sign(add(x,(-mid,0)),T,D)>=0:lo=mid
        else:hi=mid
    return lo
def shifts(x,T,D):
    n=floor(x,T,D)
    return [add(x,(-n,0))] if x[1] else [(0,0),(1,0)]
def binary(x):return x[1]==0 and x[0] in (0,1)

def main():
    if not __debug__:raise RuntimeError('no -O')
    start=time.monotonic();triples=[]
    # Independent rectangular T/N enumeration, not the searcher's D loop.
    for T in range(-16,17):
        for N in range(-71,72):
            D=T*T-4*N
            if D>0 and math.isqrt(D)**2!=D and 144*D<(203-12*abs(T))**2:triples.append((T,N,D))
    assert len(triples)==669
    counts={'C4':0,'C5':0,'C6':0};states=[];per=[]
    for T,N,D in sorted(triples):
        a=(0,1);b=(T,-1);g3=-floor(a,T,D);q3=-floor(b,T,D);A3=add(a,(g3,0));B3=add(b,(q3,0));assert add(A3,B3)==(1,0)
        for c in range(-6,7):
            before=counts.copy()
            for u,v in [(0,0),(1,0),(0,1)]:
                counts['C4']+=1;g2=u-c;q2=v-c
                for A5 in shifts(scale(a,g2),T,D):
                    for B5 in shifts(scale(b,q2),T,D):
                        if not binary(add(A5,B5,scale(B3,u),scale(A3,v))):continue
                        counts['C5']+=1
                        for A6 in shifts(scale(a,g3),T,D):
                            for B6 in shifts(scale(b,q3),T,D):
                                if not binary(add(A6,B6,mul(A3,B3,T,N))):continue
                                counts['C6']+=1
                                G=[1,0,g2,g3,-c*g2];Q=[1,0,q2,q3,-c*q2]
                                G.append(add(A5,scale(a,-g2),(-c*g3,0))[0]);Q.append(add(B5,scale(b,-q2),(-c*q3,0))[0])
                                G.append(add(A6,scale(a,-g3),(-c*G[4],0))[0]);Q.append(add(B6,scale(b,-q3),(-c*Q[4],0))[0])
                                states.append([T,N,D,c,u,v,*G,*Q])
            per.append([T,N,D,c,*[counts[k]-before[k] for k in ('C4','C5','C6')]])
    rows=sorted(states);out={'status':'complete_all_integer_endpoint_branches','TN':669,'counts':counts,'states':rows,'per_TN_c':per,'elapsed_seconds':time.monotonic()-start,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    p=Path(__file__).with_name('s3-prefix6-independent.json');p.write_text(json.dumps(out,separators=(',',':'))+'\n');print(json.dumps({k:out[k] for k in ['status','TN','counts','elapsed_seconds']}))

if __name__=='__main__':main()
