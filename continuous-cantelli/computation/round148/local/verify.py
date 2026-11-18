#!/usr/bin/env python3
"""mpmath whole-integral and phase representation checks for Metal extrema."""
import json
import re
import statistics
import struct
import sys
from functools import lru_cache
from pathlib import Path

import mpmath as mp

mp.mp.dps=65
NORM=1/mp.sqrt(2*mp.pi)
B1=mp.mpf(1)/4
B2=-mp.exp(-mp.mpf(1)/4)/4
S=13*mp.exp(-mp.mpf(1)/2)/64
C=169*mp.exp(-mp.mpf(1)/2)/2048
MASS=mp.mpf(3)/4


def f32(x):
    return struct.unpack("=f",struct.pack("=f",x))[0]


@lru_cache(maxsize=None)
def true_integral(n,turn):
    eps=1/(mp.pi*(n+mp.mpf(turn)))
    angle=2*mp.pi*mp.mpf(turn)
    wcut=[-25,-12,-6,-3,0,3,6,12,25]
    values=[]
    for d in (1,2):
        def fn(w):
            fac=1+eps*w
            return mp.sqrt(d)*NORM*fac**(-mp.mpf("2.5"))*mp.exp(-d*w*w/(2*fac))
        def oscillation(w):
            phase=2*w/(1+mp.sqrt(1+eps*w))
            return fn(w)*mp.exp(mp.j*phase)
        values.append(mp.quad(oscillation,wcut))
    combo=B1*values[0]+B2*values[1]
    # The untruncated mass transfer is 1 + 1/(d*t). The omitted region is
    # exponentially small relative to working precision at all retained t.
    residual=eps**2*MASS + mp.re(mp.exp(mp.j*angle)*combo)
    second_vector=(combo+mp.j*S*eps)/eps**2
    second=(residual-eps*S*mp.sin(angle))/eps**2
    return eps, residual, second, second_vector


def ideal_node_error(order):
    x, w = mp.gauss_quadrature(order,"legendre")
    sums=[]
    for d in (1,2):
        got=mp.mpf(0)
        for k in range(order):
            r=mp.mpf("5.5")*(1+x[k])
            got+=11*w[k]*mp.sqrt(d)*NORM*mp.exp(-d*r*r/2)*mp.cos(r)
        sums.append(got)
    return str(B1*sums[0]+B2*sums[1])


def read_samples(data):
    out=[]
    for band in data["strata"]:
        for modes in band:
            for item,label in zip(modes,data["metrics"]):
                a=item["best"]
                if item["bad"] or not a:
                    raise ValueError("Nonfinite GPU row.")
                n=a["n"]
                turn=f32(a["turns"])
                eps,resid,second,vec=true_integral(n,turn)
                gpu_eps=f32(a["eps"])
                gpu_angle=f32(a["rawAngle"]) if "rawAngle" in a else f32(2/gpu_eps)
                sin_double=float(mp.sin(mp.mpf(gpu_angle)))
                target=float(mp.sin(2*mp.pi*mp.mpf(turn)))
                out.append(dict(band=a["band"],mode=a["mode"],metric=label,n=n,turns=turn,
                                t=float(1/eps**2),eps_gpu=gpu_eps,eps_cpu=float(eps),
                                theoretical_first=float(S*target),cpu_scaled=float(resid/eps),
                                gpu_scaled=f32(a["stableScaled"]),cpu_second=float(second),
                                gpu_second=f32(a["second"]),cpu_second_real=float(mp.re(vec)),
                                cpu_second_imag=float(mp.im(vec)),gpu_second_real=f32(a["secondReal"]),
                                gpu_second_imag=f32(a["secondImag"]),raw_scaled=f32(a["rawScaled"]),
                                raw_sin_gpu=f32(a["rawSine"]),raw_sin_reference=sin_double,
                                intended_sin=target,quad_floor_gpu=f32(a["quadFloor"]),
                                raw_angle=gpu_angle,raw_angle_recorded="rawAngle" in a))
    return out


def main():
    data=json.loads(Path(sys.argv[1]).read_text())
    rows=read_samples(data)
    result=dict(samples=rows,mpmath=mp.__version__,digits=mp.mp.dps,
                ideal_quadrature_floor={str(k):ideal_node_error(k) for k in (32,64,72)},
                max_second_error=max(abs(x["gpu_second"]-x["cpu_second"]) for x in rows),
                max_phase_only_scaled_error=max(abs(S*(x["raw_sin_gpu"]-x["intended_sin"])) for x in rows),
                max_sin_disagreement_on_float_argument=max(
                    abs(x["raw_sin_reference"]-x["raw_sin_gpu"]) for x in rows),
                angle_recorded=all(x["raw_angle_recorded"] for x in rows),
                max_float_angle_representation_error=max(
                    abs(x["raw_sin_reference"]-x["intended_sin"]) for x in rows),
                max_gpu_quadrature_floor=max(abs(x["quad_floor_gpu"]) for x in rows))
    if len(sys.argv)>3:
        used=[int(m.group(1)) for line in Path(sys.argv[3]).read_text().splitlines()
              if line.startswith("during ") and (m:=re.search(r'"Device Utilization %"=(\d+)',line))]
        result["utilization"]=dict(count=len(used),median=statistics.median(used),max=max(used))
    for k,v in result.items():
        if k!="samples":
            print(k,v,flush=True)
    for mode in range(4):
        row=next(x for x in rows if x["band"]==3 and x["mode"]==mode)
        print("cardinal",mode,"t",row["t"],"first",row["cpu_scaled"],"second",row["cpu_second"],
              "raw sine",row["raw_sin_gpu"],flush=True)
    result["max_phase_only_scaled_error"]=float(result["max_phase_only_scaled_error"])
    Path(sys.argv[2]).write_text(json.dumps(result,indent=2)+"\n")


if __name__=="__main__":
    main()
