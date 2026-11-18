#!/usr/bin/env python3
"""Independent CPU spatial and positive-Fourier checks for Metal gap couplings."""
import json
import math
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np
import scipy
from scipy.special import roots_legendre

SQRT2PI=math.sqrt(2*math.pi)


@lru_cache(maxsize=None)
def points(n):
    v=np.zeros(1,dtype=np.int64)
    for _ in range(n):
        v=np.column_stack((3*v,3*v+2)).ravel()
    return (v+.5)/3**n


def densities(x):
    return np.exp(-.5*x*(4+x))


def spatial(aDegree,xDegree,m,sigma,use_gap=False,order=24):
    a,x=points(aDegree),points(xDegree)
    omega=2*np.pi*3**m
    if use_gap:
        z,w=roots_legendre(order)
        samples=(z/(2*3**(xDegree+1)),w/2)
    else:
        samples=(np.zeros(1),np.ones(1))
    deltas,weights=samples
    spread=max(abs(deltas),default=0.)
    ga=densities(a)
    ans=0j
    for start in range(0,len(a),128):
        A=a[start:start+128]
        l=x.searchsorted(A-12*sigma-spread,side="left")
        r=x.searchsorted(A+12*sigma+spread,side="right")
        if max(r-l)==0:
            continue
        xind=l[:,None]+np.arange(max(r-l))[None,:]
        mask=xind<r[:,None]
        X=x[np.minimum(xind,len(x)-1)]
        diff=A[:,None]-X
        subtotal=0j
        for du,weight in zip(deltas,weights):
            z=diff-du
            subtotal+=weight*np.sum(mask*densities(X+du)*ga[start:start+128,None]
                                    *np.exp(-.5*(z/sigma)**2-1j*omega*z))
        ans+=subtotal
    return ans/len(a)/len(x)


def fourier(n,m,sigma):
    # Integrate the explicitly positive Fourier integrand in standard
    # Gaussian Z coordinates. Used for moderate frequencies only.
    assert m<=3
    p=points(n)
    w=densities(p)
    omega=2*np.pi*3**m
    order=12
    edges=np.linspace(-11,11,1+math.ceil(22/max(.5*sigma,1e-4)))
    z,v=roots_legendre(order)
    node=((edges[:-1,None]+edges[1:,None])/2+np.diff(edges)[:,None]*z/2).ravel()
    gw=(np.diff(edges)[:,None]*v/2).ravel()*np.exp(-node**2/2)/SQRT2PI
    ans=0.
    for j in range(0,len(node),128):
        freq=omega+node[j:j+128]/sigma
        coeff=np.exp(-1j*freq[:,None]*(p-.5))@w/len(p)
        ans+=float(np.dot(gw[j:j+128],abs(coeff)**2))
    return ans


def coupling(c):
    s,j,m,n=c["sigma"],c["j"],c["m"],c["n"]
    E,H,G=map(lambda k:complex(*c[k]),("energy","midpoint","gap"))
    omega=2*np.pi*3**m
    first=2+np.exp(-.5)/s+omega
    second=3+4*np.exp(-.5)/s+1/s**2+2*omega*(2+np.exp(-.5)/s)+omega**2
    ell=3.**(-j)
    midpoint_truncation=3.**(-2*n)*second/16
    lipschitz=ell*first/3
    second_gap=29*ell**2*second/432
    assert E.real>0
    assert abs(E-G)<lipschitz and abs(E-G)<second_gap
    return dict(**c,D=s*s*omega,phase_exponent=.5*(s*omega)**2,
                gap_error=abs(E-G),midpoint_error=abs(E-H),
                relative_gap_error=abs(E-G)/E.real,lipschitz=lipschitz,
                quadratic=second_gap,a_measure_midpoint_truncation=midpoint_truncation,
                lip_ratio=abs(E-G)/lipschitz,quadratic_ratio=abs(E-G)/second_gap,
                positivity_from_lipschitz=(E.real-2*midpoint_truncation>lipschitz).item(),
                positivity_from_quadratic=(E.real-2*midpoint_truncation>second_gap).item())


def main():
    data=json.loads(Path(sys.argv[1]).read_text())
    rows=[coupling(r) for r in data["cases"]]
    high=max(rows,key=lambda r:r["relative_gap_error"])
    sign_difference=next((r for r in rows if r["gap"][0]>0 and r["midpoint"][0]<0 and r["m"]>=3),None)
    ids={(high["m"],high["j"],high["kappa"])}
    if sign_difference:
        ids.add((sign_difference["m"],sign_difference["j"],sign_difference["kappa"]))
    if "beta" in high:
        ids.update((r["m"],r["j"],r["kappa"]) for r in rows)
    ids.update((7,j,1.4) for j in range(6,12))
    ids.update([(0,0,3.5),(0,4,.20),(1,1,1.4),(2,5,.42),(3,6,.8),(9,9,3.5),(9,12,1.4)])
    checks=[]
    for c in rows:
        m,j,n,k,s=(c[key] for key in ("m","j","n","kappa","sigma"))
        if (m,j,k) not in ids:
            continue
        a=spatial(n,n,m,s)
        h=spatial(n,j,m,s)
        g=spatial(n,j,m,s,True)
        e_plus=spatial(n+1,n+1,m,s)
        g_plus=spatial(n+1,j,m,s,True)
        delta=[abs(complex(*c[key])-z) for key,z in zip(("energy","midpoint","gap"),(a,h,g))]
        assert abs(a.imag)<1e-10
        row=dict(m=m,j=j,n=n,kappa=k,cpu_energy=[a.real,a.imag],cpu_midpoint=[h.real,h.imag],
                 cpu_gap=[g.real,g.imag],cpu_gpu_abs=delta,cpu_next_depth_energy=e_plus.real,
                 cpu_next_depth_gap=[g_plus.real,g_plus.imag])
        if m<=3:
            en=fourier(n,m,s)
            row["positive_fourier_energy"]=en
            row["positive_fourier_abs_difference"]=abs(en-a.real)
        checks.append(row)
    stats=dict(cases=len(rows),gap_negative=sum(r["gap"][0]<0 for r in rows),
               midpoint_gap_sign_disagreement=sum(r["gap"][0]*r["midpoint"][0]<0 for r in rows),
               conservative_lip_positive_count=sum(r["positivity_from_lipschitz"] for r in rows),
               conservative_quadratic_positive_count=sum(r["positivity_from_quadratic"] for r in rows),
               max_error_over_lipschitz=max(r["lip_ratio"] for r in rows),
               max_error_over_quadratic_bound=max(r["quadratic_ratio"] for r in rows),
               max_cpu_gpu_error=max(max(c["cpu_gpu_abs"]) for c in checks),
               max_spatial_vs_fourier=max(c.get("positive_fourier_abs_difference",0) for c in checks))
    by_delta={}
    for k in sorted({r["j"]-r["m"] for r in rows}):
        z=[r for r in rows if r["j"]-r["m"]==k]
        by_delta[k]=dict(cases=len(z),gap_negative=sum(r["gap"][0]<0 for r in z),
                         worst_rel_error=max(r["relative_gap_error"] for r in z),
                         median_rel_error=float(np.median([r["relative_gap_error"] for r in z])),
                         pass_lip=sum(r["positivity_from_lipschitz"] for r in z),
                         pass_quadratic=sum(r["positivity_from_quadratic"] for r in z))
    out=dict(stats=stats,by_delta=by_delta,max_relative_error=high,sign_difference=sign_difference,
              checked=checks,all_cases=rows,numpy=np.__version__,scipy=scipy.__version__)
    Path(sys.argv[2]).write_text(json.dumps(out,indent=2)+"\n")
    print(stats)
    print(by_delta)
    for check in checks:
        print(check)


if __name__=="__main__":
    main()
