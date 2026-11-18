"""High-order Hermite leakage jets with physical/log-parameter crosschecks."""
import json, math, sys, time

import mpmath as mp
from scipy.stats import qmc
import flat_core as base


def side_jet(a,p,q,ell,sign,k,log_input):
    R=mp.mpf(base.R)
    const=2*mp.log(q*R)
    q0=1/mp.sqrt(2*mp.pi)
    def f(u):
        if u<=0:
            return mp.mpf(0)
        w=const-ell*((R/u)**p-1)
        log_z=mp.log(u)-w/2
        if log_z>10:
            return mp.mpf(0)
        z=mp.exp(log_z)
        return q0*sign**k*mp.exp(-w*(k+1)/2-z*z/2-sign*a*u-u*u/2)*base.hermite(k,z)
    if log_input:
        return mp.quad(lambda t: (u:=R*mp.exp(-t))*f(u),[0,.125,.25,.5,.75,1,1.5,2,3,4])
    return mp.quad(f,[0,R/32,R/16,R/8,R/4,R/2,3*R/4,7*R/8,R])


def worker(rows,smoke):
    mp.mp.dps=57
    start=base.utc()
    cpu=time.process_time()
    cases=[]
    for index,a,p,q,ell,qq,ll in rows:
        args=[(1,q,ell)]+([(-1,qq,ll)] if qq else [])
        a,p=mp.mpf(a),mp.mpf(p)
        jets=[]
        for k in (0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64):
            phys,logscale=[],[]
            for sign,q_j,l_j in args:
                q_j,l_j=mp.mpf(q_j),mp.mpf(l_j)
                phys.append(side_jet(a,p,q_j,l_j,sign,k,False))
                logscale.append(side_jet(a,p,q_j,l_j,sign,k,True))
            denom=max(1,*map(abs,phys),*map(abs,logscale))
            cross=max(abs(x-y) for x,y in zip(phys,logscale))/denom
            if cross>mp.mpf("1e-20"):
                raise ArithmeticError((index,k,mp.nstr(cross,20)))
            coefficient=mp.fsum(phys)/mp.factorial(k)
            jets.append({"k":k,"taylor_coefficient":mp.nstr(coefficient,30),
                         "root_growth":mp.nstr(abs(coefficient)**(mp.mpf(1)/k),22) if k else None,
                         "mp_cross_scaled":mp.nstr(cross,20)})
        cases.append({"i":index,"a":float(a),"p":float(p),"q_plus":q,"ell_plus":ell,
                      "q_minus":qq,"ell_minus":ll,"jets":jets})
    return {"pid":__import__("os").getpid(),"start":start,"finish":base.utc(),
            "cpu_s":time.process_time()-cpu,"cases":cases}


if __name__=="__main__":
    base.DPS=57
    base.SEED=301109
    base.worker=worker
    base.main()
