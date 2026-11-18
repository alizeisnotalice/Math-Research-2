#!/usr/bin/env python3
"""Third S13 check using the parity of C4 rather than intersecting d shifts."""
from argparse import ArgumentParser
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from math import isqrt
from pathlib import Path
import json
import hashlib
import sys
import time

Z=(0,0)
I=(1,0)
BOUNDS={2:(17,1),3:(65,2),4:(257,4),5:(1025,8),6:(4097,16),7:(16385,32),8:(65537,64)}
RAD={2:41,4:120,5:92,6:41,7:9}


def add(a,b):
    return a[0]+b[0],a[1]+b[1]


def scale(s,a):
    return s*a[0],s*a[1]


def ceil(n,d):
    return -((-n)//d)


class Field:
    def __init__(self, D):
        self.D=D
        self.t=int(D%4==1)
        self.tau=1 if self.t else 2
        self.norm_carry=(D-1)//4 if self.t else D

    def mul(self,A,B):
        x,y=A
        u,v=B
        return (x*u+y*v*self.norm_carry,x*v+y*u+y*v*self.t)

    def conj(self,A):
        x,y=A
        return x+self.t*y,-y

    def sign(self,A):
        u,v=A
        base=2*u+self.t*v
        d=self.tau*v
        if base>=0 and d>=0:return int(base+d>0)
        if base<=0 and d<=0:return -int(base+d<0)
        k=base*base-d*d*self.D
        assert k
        return (2*int(k>0)-1)*(2*int(base>0)-1)

    def house(self,A,k):
        p,q=BOUNDS[k]
        u,v=A
        R=2*p-q*abs(2*u+self.t*v)
        return R>0 and R*R>self.D*(q*self.tau*v)**2

    def integer_candidates(self,base,k):
        num,den=BOUNDS[k]
        u,v=base
        rad=(2*num-isqrt(self.D*(den*self.tau*v)**2)-1)//den
        if rad<0:return range(0)
        t=2*u+self.t*v
        return range(max(-RAD[k],ceil(t-rad,2*k)),min(RAD[k],(t+rad)//(2*k))+1)

    def around(self,num):
        """All nonzero-radical coefficients with symmetric strict house num/2."""
        D=self.D
        for s in range(1,isqrt((num*num-1)//(self.tau*self.tau*D))+1):
            gap=num-isqrt((self.tau*s)**2*D)-1
            for m in range(ceil(-gap-self.t*s,2),(gap-self.t*s)//2+1):
                yield m,s

    def c_candidates(self,base):
        """Exact two-embedding constraints on c from house and s3."""
        u,v=base
        D=self.D
        U,V=2*u+self.t*v,self.tau*v
        r=isqrt((65**2-1)//D)
        smax=isqrt((185**2-1)//(D*self.tau**2))
        for s in range(max(-smax,ceil(V-r,3*self.tau)),min(smax,(V+r)//(3*self.tau))+1):
            if not s:continue
            radical=abs(V-3*self.tau*s)
            newton_radius=64-isqrt(radical*radical*D)
            Tbound=184-isqrt((self.tau*s)**2*D)
            low=max(-Tbound,ceil(U-newton_radius,3))
            top=min(Tbound,(U+newton_radius)//3)
            for m in range(ceil(low-self.t*s,2),(top-self.t*s)//2+1):
                yield m,s

    def quotient(self,base,zero_only):
        m,n=base
        if zero_only:
            return ((-m,Z),) if not n else ()
        if not n:
            return ((-m,Z),(1-m,I))
        u,v=2*m+self.t*n,self.tau*n
        w=isqrt(v*v*self.D)
        if v<0:w=-w-1
        h=-((u+w)//2)
        out=(m+h,n)
        assert self.sign(out)>0 and self.sign(add(I,scale(-1,out)))>0
        return ((h,out),)

    def conv(self,p,q,k):
        result=Z
        for i in range(max(0,k-len(q)+1),min(len(p)-1,k)+1):
            result=add(result,self.mul(p[i],q[k-i]))
        return result

    def power_base(self,p,s,k):
        return scale(-1,self.conv([Z,*p[1:k]],s,k))

    def step(self,p,states,k,force_left=False,force_right=False):
        pb=tuple(map(self.conj,p))
        answer=[]
        for G,Q,A,B in states:
            ua=self.conv(p,[(g,0) for g in G],k)
            ub=self.conv(pb,[(q,0) for q in Q],k)
            cross=self.conv(A,B,k)
            for g,an in self.quotient(ua,force_left):
                for q,bn in self.quotient(ub,force_right):
                    if add(cross,add(an,bn)) in (Z,I):
                        answer.append((G+(g,),Q+(q,),A+(an,),B+(bn,)))
        return answer


INITIAL=[((1,),(1,),(I,),(I,))]


def field(D):
    F=Field(D)
    counts=Counter()
    records=[]
    for a in F.around(20):
        counts['a_cores']+=1
        pa=[I,a]
        a1=F.step(pa,INITIAL,1)
        row2=F.step(pa+[Z],a1,2)
        if not row2:continue
        counts['a_row2']+=1
        assert len(a1)==len(row2)==1
        a_sq=F.mul(a,a)
        s1=scale(-1,a)
        G,Q,A,B=row2[0]
        for b in F.integer_candidates(a_sq,2):
            counts['newton2_c']+=1
            P2=[*pa,(b,0)]
            s2=add(a_sq,(-2*b,0))
            U2=[(G[:2]+(G[2]-b,),Q[:2]+(Q[2]-b,),A,B)]
            powers=[(8,0),s1,s2]
            base3=F.power_base(P2,powers,3)
            for c in F.c_candidates(base3):
                counts['newton3_core']+=1
                r,s=a[1],c[1]
                g2,q2=G[2]-b,Q[2]-b
                if s*g2%r or s*q2%r:continue
                counts['row5_projection_divisibility']+=1
                P3=[*P2,c]
                Pb3=list(map(F.conj,P3))
                row3=F.step(P3,U2,3)
                counts['row3']+=len(row3)
                counts['row4_normalized']+=len(F.step([*P3,Z],row3,4))
                for W,V,A3,B3 in row3:
                    g4,q4=-s*g2//r,-s*q2//r
                    ua=F.conv(P3,[(g,0) for g in W+(g4,)],4)
                    ub=F.conv(Pb3,[(q,0) for q in V+(q4,)],4)
                    cross=F.conv(A3,B3,4)
                    assert F.sign(cross)>0 and F.sign(add(I,scale(-1,cross)))>0
                    constant=add(add(ua,ub),cross)
                    assert constant[1]==0
                    d=(2-constant[0])//2
                    ax,bx=add(ua,(d,0)),add(ub,(d,0))
                    if not all(F.sign(y)>=0 and F.sign(add(I,scale(-1,y)))>=0 for y in (ax,bx)):continue
                    C4=constant[0]+2*d
                    assert C4 in (1,2)
                    counts['shared_unit_d']+=1
                    counts[f'C4_value_{C4}']+=1
                    key=(D,2*a[0]+F.t*r,r,b,2*c[0]+F.t*s,s,d)
                    records.append((key,C4,ax,bx))
    return counts,records


def main():
    if not __debug__:raise SystemExit('requires assertions')
    p=ArgumentParser(description=__doc__)
    p.add_argument('--reference',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--workers',type=int,default=4)
    a=p.parse_args()
    start=time.monotonic()
    ds=[D for D in range(2,400) if all(D%(k*k)!=0 for k in range(2,isqrt(D)+1))]
    result=Counter()
    rows=[]
    with ProcessPoolExecutor(a.workers) as ex:
        for counts,local in ex.map(field,ds):
            result.update(counts)
            rows.extend(local)
    result=Counter({name:num for name,num in result.items() if num>0})
    payload=b''.join(json.dumps(row,separators=(',',':')).encode()+b'\n' for row in sorted(rows))
    original=a.reference.read_bytes()
    assert payload==original
    assert len(rows)==866==result['shared_unit_d']==result['C4_value_2'] and 'C4_value_1' not in result
    with a.out.with_suffix('.rows.jsonl').open('xb') as f:f.write(payload)
    receipt={'status':'verified','counts':dict(result),'byte_identical':True,'output_sha256':hashlib.sha256(payload).hexdigest(),
             'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             'python':sys.version,'workers':a.workers,'seconds':time.monotonic()-start}
    with a.out.open('x') as f:f.write(json.dumps(receipt,sort_keys=True,indent=2)+'\n')
    print(json.dumps(receipt,sort_keys=True))


if __name__=='__main__':main()
