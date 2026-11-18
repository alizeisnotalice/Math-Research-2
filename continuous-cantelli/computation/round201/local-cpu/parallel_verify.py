import argparse
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import platform
import subprocess
from time import monotonic

import mpmath as m

ROOT=Path(__file__).parent


def utc():
    return datetime.now(timezone.utc).isoformat()


def c(k):
    return m.power(2,k-m.mpf('.5'))*m.gamma(k+m.mpf('.5'))/m.pi


def factor(k):
    return m.gamma(1-k)/(k*m.power(2,k))


def kernel(k,y):
    z=y*y/2
    return c(k)*m.exp(-z)*m.hyp1f1(-k,m.mpf('.5'),z)


def marchaud(k,y):
    z=y*y/2
    phi=m.exp(-z)/m.sqrt(2*m.pi)
    a=1-k
    b=k+m.mpf('.5')

    def early(u):
        if u==0:
            return (m.mpf('.5')-z)*phi/a
        t=u**(1/a)
        return -m.expm1(-m.log1p(t)/2+z*t/(1+t))/t*phi/a

    def late(u):
        v=u**(1/b)
        return m.exp(-z*v/(1+v))/m.sqrt(2*m.pi*(1+v))/b

    split=m.power(1/(1+y*y),b)
    cuts=sorted({m.mpf(0),split/4,split,min(m.mpf(1),4*split),m.mpf(1)})
    return m.quad(early,[0,.5,1])+phi/k-m.quad(late,cuts)


def fourier(k,y):
    q=lambda s:s**(2*k)*m.exp(-s*s/2)*m.cos(y*s)/m.pi
    return m.quad(q,[0,1,2,3,5,7,10,14,20])


def normalized_u(k,z):
    """U_k times exp(-z) z^(k+1/2) Gamma(1-k)/sqrt(pi)."""
    b=k+m.mpf('.5')
    a=1-k

    def right(v):
        u=v**(1/b)
        return (m.exp(-u)-m.exp(-z))*(1-u/z)**(-1-k)/b

    def left(v):
        if v==0:
            return z/a
        t=v**(1/a)
        return m.expm1(z*t)/t*(1-t)**(b-1)/a

    raw=[m.mpf(q) for q in [0,1,4,12,30,60,120] if q<z/2]+[z/2]
    cuts=[u**b for u in raw]
    return (m.quad(right,cuts)+z**b*m.exp(-z)*m.quad(left,[0,(m.mpf('.5')**a)/2,m.mpf('.5')**a]))/m.gamma(b)


def root_row(l):
    l=m.mpf(l)
    k=m.exp(-l)
    ell=m.log(l/m.pi)
    lead=l+ell/2
    second=lead+(ell/4-m.mpf('.5'))/l
    third=second+(-ell**2/16+3*ell/8-m.mpf('7/8'))/(l*l)

    def log_equation(z):
        return z-l+m.log(m.pi)/2-m.loggamma(1-k)-(k+m.mpf('.5'))*m.log(z)+m.log(normalized_u(k,z))

    root=m.findroot(log_equation,(lead-m.mpf('.12'),lead+m.mpf('.12')),tol=m.mpf('1e-39'))
    assert abs(log_equation(root))<m.mpf('2e-37')
    detail=''
    if l<=160:
        with m.workdps(65):
            kk=m.exp(-l)
            a,b=root-m.mpf('.03'),root+m.mpf('.03')
            assert m.hyp1f1(-kk,m.mpf('.5'),a)>0>m.hyp1f1(-kk,m.mpf('.5'),b)
            for _ in range(127):
                mid=(a+b)/2
                if m.hyp1f1(-kk,m.mpf('.5'),mid)>0:
                    a=mid
                else:
                    b=mid
            gap=abs(root-(a+b)/2)
        assert gap<m.mpf('1e-34')
        detail=m.nstr(gap,9)
    lambert=-m.lambertw(-2*m.pi*k*k,-1)/2 if l>3 else m.nan
    return dict(L=str(l),kappa_mantissa=m.nstr(k,11),z=m.nstr(root,27),
                scaled_leading_error=m.nstr(l*(root-lead),17),
                scaled_second_error=m.nstr(l*(root-second),17),
                scaled_third_error=m.nstr(l*l*(root-third),17),
                comparison_hypergeometric=detail,lambert_difference=m.nstr(root-lambert,12))


def worker(job):
    shard,count=job
    m.mp.dps=51
    stamped=utc()
    t=monotonic()
    direct=[]
    tails=[]
    for j in range(count):
        i=12*j+shard
        k=m.mpf(12+(i*71+23)%960)/1000
        y=m.mpf((i*127+13)%1100)/100
        kn=kernel(k,y)
        ma=marchaud(k,y)
        delta=abs(kn*factor(k)-ma)
        allowance=m.mpf('2e-39')*(1+abs(ma))
        assert delta<allowance,(k,y,delta)
        err=None
        if j%5==0:
            err=abs(kn-fourier(k,y))
            assert err<m.mpf('1e-40'),(k,y,err)
        direct.append(dict(kappa=str(k),y=str(y),kernel=m.nstr(kn,25),
                           marchaud_delta=m.nstr(delta,12),fourier_delta=m.nstr(err,12) if err is not None else None))
        if j%3==0:
            y=m.mpf(15+11*j+shard)
            ma=marchaud(k,y)
            cn=2**k*m.gamma(k+m.mpf('.5'))/m.sqrt(m.pi)
            ratio=-ma*y**(1+2*k)/cn
            next_coeff=(2*k+1)*(2*k+2)/2
            second_coeff=(2*k+1)*(2*k+2)*(2*k+3)*(2*k+4)/8
            error=(ratio-1-next_coeff/y**2)
            tails.append(dict(kappa=str(k),y=str(y),ratio=m.nstr(ratio,15),
                        y4_remainder=m.nstr(y**4*error,13),next_coefficient=m.nstr(second_coeff,13)))
            assert not ma>=0
    Ls=['3','4','6','9','12','16','20','28','40','70','100','160',
        '250','400','700','1000','1600','2500','5000','10000','25000','50000','100000','1000000']
    roots=[root_row(l) for l in Ls[shard::12]]
    result=dict(pid=os.getpid(),shard=shard,started=stamped,finished=utc(),seconds=monotonic()-t,
                direct=direct,tails=tails,roots=roots)
    (ROOT/f'shard-{shard:02d}.json').write_text(json.dumps(result,indent=2)+'\n')
    return {k:v for k,v in result.items() if k not in ('direct','tails','roots')}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--count',type=int,default=14)
    cfg=ap.parse_args()
    launch=utc()
    with ProcessPoolExecutor(max_workers=12) as pool:
        jobs=[pool.submit(worker,(k,cfg.count)) for k in range(12)]
        procs=pool._processes
        pids=list(procs)
        print(json.dumps(dict(utc=launch,parent=os.getpid(),worker_pids=pids)),flush=True)
        outcomes=[j.result() for j in jobs]
    alive=subprocess.run(['ps','-p',','.join(map(str,pids)),'-o','pid='],capture_output=True,text=True)
    rows=[json.loads(p.read_text()) for p in ROOT.glob('shard-*.json')]
    direct=[v for q in rows for v in q['direct']]
    tails=[v for q in rows for v in q['tails']]
    roots=sorted([v for q in rows for v in q['roots']],key=lambda r:float(r['L']))
    result=dict(start=launch,end=utc(),parent=os.getpid(),python=platform.python_version(),mpmath=m.__version__,
                workers=outcomes,n_direct=len(direct),n_tails=len(tails),n_roots=len(roots),
                n_fourier=sum(r['fourier_delta'] is not None for r in direct),
                maximum_marchaud_difference=max(float(r['marchaud_delta']) for r in direct),
                maximum_fourier_difference=max(float(r['fourier_delta']) for r in direct if r['fourier_delta'] is not None),
                maximum_root_formula_difference=max(float(r['comparison_hypergeometric'])
                                                    for r in roots if r['comparison_hypergeometric']),
                roots=roots,exit_codes={pid:proc.exitcode for pid,proc in procs.items()},remaining=alive.stdout.splitlines())
    (ROOT/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['workers','roots']}),flush=True)


if __name__=='__main__':
    main()
