"""A certified lower integral for a uniform radial entropy dual."""
import argparse,json,os,platform,time
from decimal import Decimal as D
from fractions import Fraction as F
from pathlib import Path
from interval import Box,NEAR,LO,HI,ZERO,ONE,TWO


def atan_bounds(x,n=12):
    s=sum((-1)**i*x**(2*i+1)/(2*i+1) for i in range(n))
    return s,s+x**(2*n+1)/(2*n+1)


a,b=atan_bounds(F(1,5)),atan_bounds(F(1,239))
pi_lo,pi_hi=16*a[0]-4*b[1],16*a[1]-4*b[0]
assert F(314159,100000)<pi_lo<pi_hi<F(314160,100000)
C=Box('1.2533','1.2534')
assert 2*F(12533,10000)**2<pi_lo and pi_hi<2*F(6267,5000)**2


def log(x):
    assert x.a>0
    return Box(NEAR.next_minus(NEAR.ln(x.a)),NEAR.next_plus(NEAR.ln(x.b)))


def hyper(x):
    e=x.exp();en=(ZERO-x).exp()
    return (e-en)/TWO,(e+en)/TWO


def kernel(z,x):
    Z,T=Box(z),Box(x)
    U=Z*T
    sh,ch=hyper(U)
    return C*sh/Z, C*(U*ch-sh)/Z.pow(3)


def grid(step):
    nodes=[D('.02')+i*step for i in range(1+int(D('11.98')/step))]
    if nodes[-1]<D(12):nodes.append(D(12))
    assert nodes[0]==D('.02') and nodes[-1]==12
    logcosh=[log(hyper(Box('1.8')*Box(x))[1]) for x in nodes]
    gamma=[(ZERO-Box(x).pow(2)/TWO).exp()/C for x in nodes]
    return nodes,logcosh,gamma


def integrate(zlo,zhi,data):
    q=Box(zlo,zhi).pow(2)/TWO
    scale=Box('1.62')-q
    offset=ONE+log(ONE-q)
    t=Box('1.5')
    nodes,lc,gauss=data
    lows=[kernel(zlo,x) for x in nodes]
    highs=[kernel(zhi,x) for x in nodes]
    total=D(0)
    for i in range(len(nodes)-1):
        a,b=nodes[i:i+2]
        k=Box(lows[i][0].a,highs[i+1][0].b)
        kq=Box(lows[i][1].a,highs[i+1][1].b)
        H=k*(Box(lc[i].a,lc[i+1].b)-log(k)+offset)-scale*kq
        h=H.a
        m=h if h<=t.a else (t*(ONE+log(Box(h)/t))).a
        density=gauss[i+1].a if m>=0 else gauss[i].b
        total=LO.add(total,LO.multiply(LO.multiply(m,density),b-a))
    # Tail bound: integral from 12 of x^3 exp(-.45*x*x), validated analytically.
    tail=(Box('65.8')*(ZERO-Box('64.8')).exp()/Box('.405')).b
    exact=LO.subtract(LO.subtract(total,tail),(ONE+q/TWO).b)
    return dict(z=[str(zlo),str(zhi)],panels=len(nodes)-1,tail_upper=str(tail),margin_lower=str(exact))


def main():
    p=argparse.ArgumentParser();p.add_argument('shard',type=int);p.add_argument('out',type=Path)
    p.add_argument('--workers',type=int,default=15);p.add_argument('--dx',type=D,default=D('.01'))
    p.add_argument('--dz',type=D,default=D('.0025'));args=p.parse_args()
    start=time.time();cpu=time.process_time()
    results=[];dat=grid(args.dx)
    last=int(D('.59')/args.dz)
    assert D(last)*args.dz==D('.59')
    for i in range(args.shard,last,args.workers):
        zl=D('.01')+i*args.dz
        results.append(integrate(zl,zl+args.dz,dat))
    args.out.write_text(json.dumps(dict(shard=args.shard,start=start,end=time.time(),pid=os.getpid(),
        cpu=time.process_time()-cpu,python=platform.python_version(),results=results),indent=2)+'\n')
    assert all(D(x['margin_lower'])>0 for x in results)

if __name__=='__main__':main()
