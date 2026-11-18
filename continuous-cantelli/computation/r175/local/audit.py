#!/usr/bin/env python3
"""Independent adaptive SciPy audit of every retained native peak extremum."""
import json
import re
import statistics
import struct
import sys
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import quad_vec

C=1/np.sqrt(2*np.pi)


def f32(x):
    return struct.unpack("=f", struct.pack("=f",x))[0]


def measure(shape):
    c,L,s,w1,w2,h1,h2,l1,r1,l2,r2=(f32(shape[k]) for k in
        ("c","L","s","w1","w2","h1","h2","l1","r1","l2","r2"))

    def density(u,v):
        x=c+L*(u-.5)
        return C*np.exp(-.5*x*x/v)/np.sqrt(v)

    def kernel(u,h,v):
        return density(u,v)*np.array([1.,u,u*u-h,np.sin(np.pi*u)*np.exp(np.pi*np.pi*h/2)])

    def integrate(fn,l,r):
        ans,err=quad_vec(fn,l,r,epsabs=4e-13,epsrel=4e-13)
        assert err<1e-9
        return ans

    R=integrate(lambda u:kernel(u,0.,s),0.,1.)
    P=np.zeros(4)
    for start,w,h,lo,hi in ((0.,w1,h1,l1,r1),(1-w2,w2,h2,l2,r2)):
        for a,b,ha,hb in ((0.,lo,0.,h),(lo,hi,h,h),(hi,1.,h,0.)):
            P+=integrate(lambda t:w*kernel(start+w*t,ha+(hb-ha)*(t-a)/(b-a),1.),a,b)
    d0,d1,d2,sine=R-P
    return dict(D0=d0,D1=d1,D2=d2,sine=sine,parab=d1-d2,det=d0*d2-d1*d1,
                maxValue=s+L*L*max(h1,h2),width_gap=L*(1-w1-w2))


def main():
    data=json.loads(Path(sys.argv[1]).read_text())
    selected=[]
    for k,record in enumerate(data["strata"]):
        gpu=record["best"]
        CPU=measure(gpu["shape"])
        delta={x:f32(gpu[x])-CPU[x] for x in ("D0","D1","D2","parab","det","sine")}
        assert CPU["D0"]>0 and CPU["parab"]>0 and CPU["det"]>0 and CPU["sine"]<0
        assert CPU["width_gap"]>0
        selected.append(dict(mode=k,cpu=CPU,gpu_minus_cpu=delta,indices=[gpu["index"],gpu["salt"]],
                             shape=gpu["shape"],tested=record["tested"],witness=record["witness"],
                             bad=record["bad"]))
    summary=dict(samples=selected,numpy=np.__version__,scipy=scipy.__version__,
                 maximum_abs_error={k:max(abs(p["gpu_minus_cpu"][k]) for p in selected)
                                    for k in selected[0]["gpu_minus_cpu"]})
    if len(sys.argv)>3:
        used=[int(m.group(1)) for line in Path(sys.argv[3]).read_text().splitlines()
              if line.startswith("during ") and (m:=re.search(r'"Device Utilization %"=(\d+)',line))]
        summary["ioreg_during"]=dict(count=len(used),median=statistics.median(used),max=max(used))
    Path(sys.argv[2]).write_text(json.dumps(summary,indent=2)+"\n")
    print("maximum abs differences",summary["maximum_abs_error"])
    print("ioreg",summary.get("ioreg_during"))
    for q in selected:
        print(q["mode"],q["cpu"]["maxValue"],q["cpu"]["width_gap"],q["cpu"]["D0"],
              q["cpu"]["det"],q["cpu"]["parab"],q["cpu"]["sine"])


if __name__=="__main__":
    main()
