"""All-index exact rational bracketing at delta=2^-14, split across workers."""
import argparse,json,math,os,time
from fractions import Fraction as F
from pathlib import Path

Q=10**27
DELTA=F(1,2**14)
S=24
LAST=math.ceil(F(252,100)/DELTA)
C=[(DELTA/2)**s/F(math.factorial(s+1)) for s in range(S+1)]


def cap_bracket(k):
    z=1;l=F(0);u=F(0)
    m=min(k//2,S)
    for s in range(1,m+1):
        z*=(k-2*s+2)*(k-2*s+1)
        sq=math.isqrt(z*Q**2)
        l+=C[s]*F(sq,Q)
        u+=C[s]*F(sq+1,Q)
    if m<k//2:
        x=DELTA*k/2
        u+=x**(m+1)/math.factorial(m+2)/(1-x/(m+3))
    return l,u


def main():
    p=argparse.ArgumentParser();p.add_argument('shard',type=int);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    start=time.monotonic();cpu=time.process_time()
    i=LAST*a.shard//15
    j=LAST*(a.shard+1)//15
    l=u=0
    for k in range(i,j):
        x,y=cap_bracket(k)
        l+=math.floor(Q*max(F(0),1-y*y))
        u+=math.ceil(Q*max(F(0),1-x*x))
    lower,_=cap_bracket(LAST)
    assert lower>1
    a.out.write_text(json.dumps(dict(shard=a.shard,pid=os.getpid(),range=[i,j],scale=Q,lower=l,upper=u,
       terminal_index=LAST,terminal_above_one=True,elapsed=time.monotonic()-start,cpu=time.process_time()-cpu),indent=2)+"\n")
if __name__=='__main__': main()
