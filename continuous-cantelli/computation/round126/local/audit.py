#!/usr/bin/env python3
"""Independent double integrations on the retained Metal evaluation parameters."""
import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
import json
import math
import re
import statistics
import struct
import sys
from pathlib import Path

from scipy.integrate import quad
from scipy.special import logsumexp
import numpy as np

from fit import NORM, quadrature, transform

NAMES = ("real_pos","real_neg","imag","complex","tail_pos","tail_neg","Q","far_plus","far_minus")


def f32(x):
    return struct.unpack("=f", struct.pack("=f", x))[0]


def cells(grid, b):
    return zip(grid[:-1],grid[1:],b[:-1],b[1:])


def q_value(grid,b,t,y,eps):
    def p(x,v):
        return NORM/math.sqrt(v)*math.exp(-x*x/(2*v))
    acc = p(y,t+eps) if t<1 else 0.
    for l,r,v0,v1 in cells(grid,b):
        if v0==v1==0:
            continue
        if t>=1:
            acc += quad(lambda x:p(x,1)*p(y-x,t-1+eps),l,r,epsabs=1e-13)[0]
        dv = v1-v0
        l0,r0=l,r
        if dv>0:
            r0=max(l,min(r,l+(r-l)*(t-1-v0)/dv))
        elif dv<0:
            l0=max(l,min(r,l+(r-l)*(t-1-v0)/dv))
        elif v0>t-1:
            continue
        if r0>l0:
            acc -= quad(lambda x:p(x,1)*p(y-x,t-1-(v0+dv*(x-l)/(r-l))+eps),
                        l0,r0,epsabs=1e-13)[0]
    return acc


def tail_parts(grid, b, t):
    q=t*t/2
    logs=[[],[]]
    for l,r,b0,b1 in cells(grid,b):
        sl=(b1-b0)/(r-l)
        for j,side in enumerate((-1,1)):
            if max(side*b0,side*b1)<=0:
                continue
            a=l if side*b0>=0 else l-b0/sl
            c=r if side*b1>=0 else l-b0/sl
            slope=t-q*sl if side<0 else t
            v=min(c,max(slope,a))
            v0=t*v-v*v/2+(-q*(b0+sl*(v-l)) if side<0 else 0)
            for sg,length in ((-1,v-a),(1,c-v)):
                if length<=0:
                    continue
                decay=abs(slope-v)
                k=max(1.,decay)
                width=min(k*length,100.)
                def inside(z):
                    x=v+sg*z/k
                    d=x-v
                    u=max(0,q*side*(b0+sl*(x-l)))
                    return math.exp((slope-v)*d-d*d/2)*(-math.expm1(-u))
                value=quad(inside,0,width,epsabs=3e-13,epsrel=3e-12)[0]
                if value>0:
                    logs[j].append(v0+math.log(NORM*value/k))
    return [float(logsumexp(row)) if row else -1e20 for row in logs]


def U_tests(grid, b):
    x,w,basis=quadrature(grid,16)
    b=basis@b
    w=w*np.exp(-x*x/2)
    z  = (np.polynomial.legendre.leggauss(20)[0]+1)/2
    ws = np.polynomial.legendre.leggauss(20)[1]/2
    variances=1-b[:,None]*z
    diffs=np.array([0., .5,1.,2.,3.,-0.5,-1.,-2.,-3.])[:,None,None]-x[None,:,None]
    p=NORM/np.sqrt(variances)*np.exp(-.5*diffs*diffs/variances)
    vals=(p@ws)@(-w*b)
    return [float(x) for x in vals]


def main():
    data=json.loads(Path(sys.argv[1]).read_text())
    audit=json.loads(Path("candidates.json").read_text())
    rows=[]
    errors={key:0. for key in NAMES}
    for candidate,modes in zip(audit["candidates"],data["byProfile"],strict=True):
        B=np.array(candidate["offset"],dtype=np.float32).astype(float)
        grid=np.linspace(-candidate["R"],candidate["R"],len(B))
        x,w,basis=quadrature(grid,32)
        report={}
        for key,row in zip(NAMES,modes,strict=True):
            g=row["best"]; assert g is not None
            p,q=f32(g["p"]),f32(g["q"])
            a,b=f32(g["a"]),f32(g["b"])
            if key in NAMES[:4]:
                if key in NAMES[:2]:
                    assert all(p!=t for t in candidate["training"])
                val=transform([complex(p,q)],x,w,basis,B)[0]
                result=[float(val.real),float(val.imag)]
            elif key=="Q":
                result=[q_value(grid,B,p,q,f32(g["aux"])),0.]
            else:
                result=tail_parts(grid,B,p)
            delta=max(abs(result[0]-a),abs(result[1]-b))
            errors[key]=max(errors[key],delta)
            report[key]=dict(p=p,q=q,aux=f32(g["aux"]),cpu=result,gpu=[a,b],abs_error=delta,
                             tested=row["tested"],gpu_negative=row["negative"],
                             nonfinite=row["nonfinite"])
        corner=np.max(np.abs(np.diff(B)))/candidate["dx"]
        a0=B[np.flatnonzero(grid==0)[0]]
        a1=B[np.flatnonzero(grid==1)[0]]
        assert a0==-0.25 and a1==0 and corner<=candidate["common_lip"]+1e-6
        assert np.max(abs(B))<=.5 and B[0]==0 and B[-1]==0
        rows.append(dict(case=candidate["case"],layer=candidate["layer"],R=candidate["R"],
                         fitted_max=candidate["max_abs_training"],observed_Lipschitz=corner,
                         full_interval_U=U_tests(grid,B), modes=report))
    result=dict(errors=errors,results=rows,
                count=sum(report["tested"] for row in rows for report in row["modes"].values()))
    if len(sys.argv)>3:
        used=[int(m.group(1)) for line in Path(sys.argv[3]).read_text().splitlines()
              if line.startswith("during ") and (m:=re.search(r'"Device Utilization %"=(\d+)',line))]
        result["utilization"]=dict(count=len(used),median=statistics.median(used),highest=max(used))
    Path(sys.argv[2]).write_text(json.dumps(result,indent=2)+"\n")
    print("maximum discrepancies",errors)
    for r in rows:
        if r["layer"]==6:
            modes=r["modes"]
            print("final",r["case"],"fit",r["fitted_max"],
                  {k:f'{complex(*modes[k]["cpu"]):.6g}' for k in NAMES[:4]},
                  "Q min",modes["Q"]["cpu"][0],"U max",max(abs(z) for z in r["full_interval_U"]),
                  "far margins",modes["far_plus"]["cpu"][0]-modes["far_plus"]["cpu"][1],
                  modes["far_minus"]["cpu"][0]-modes["far_minus"]["cpu"][1])


if __name__=="__main__":
    main()
