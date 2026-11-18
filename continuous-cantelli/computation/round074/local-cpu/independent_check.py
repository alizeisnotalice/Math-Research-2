"""Explicit physical-gap enumeration against product-form truncated sums."""
import math
from fractions import Fraction
import json
from pathlib import Path

import mpmath as mp
from scipy.integrate import quad

from cantor_fourier import calc


def direct(b,A,alpha,c,n,depth):
    intervals=[(Fraction(0),Fraction(1))]
    total=0j
    correction=0.
    gaps=0
    for _ in range(depth):
        next_intervals=[]
        for x,y in intervals:
            step=(y-x)/b
            next_intervals.extend((x+j*step,x+(j+1)*step) for j in A)
            for d,e in zip(A,A[1:]):
                if e==d+1:
                    continue
                a,bound=float(x+(d+1)*step),float(x+e*step)
                omega=2*math.pi*n
                half=(bound-a)/2
                coeff=2*math.pi**2*n*n*c
                cap=min(110.,coeff*half**alpha)
                physical=(cap/coeff)**(1/alpha)
                power=math.ceil(6*alpha)
                jac_power=power/alpha
                def value(s):
                    offset=physical*s**jac_power
                    weight=physical*jac_power*s**(jac_power-1)*math.exp(-cap*s**power)
                    left,right=omega*(a+offset),omega*(bound-offset)
                    return (weight*(math.cos(left)+math.cos(right)),
                            weight*(math.sin(left)+math.sin(right)),2*weight)
                vals=[quad(lambda t:value(t)[j],0,1,epsabs=1e-14,limit=150)[0] for j in range(3)]
                plain=complex(math.cos(omega*bound),math.sin(omega*bound))-complex(math.cos(omega*a),math.sin(omega*a))
                total+=complex(vals[0],vals[1])-plain/(1j*omega)
                correction+=vals[2]-(bound-a)
                gaps+=1
        intervals=next_intervals
    return total,1+correction,gaps


def main():
    mp.mp.dps=60
    nodes,weights=mp.gauss_quadrature(68,"legendre")
    tests=((3,(0,2),.5,.3,13,5),(5,(0,1,4),1.2,.7,37,4),
           (7,(0,2,4,6),2.,.13,49,4),(8,(0,1,3,7),3.1,.66,63,4))
    out=[]
    for b,A,alpha,c,n,depth in tests:
        exact=calc(b,A,alpha,c,n,True,(tuple(nodes),tuple(weights)),fixed_depth=depth)
        h,m,gaps=direct(b,A,alpha,c,n,depth)
        diff=max(abs(complex(exact[0])-h),abs(float(exact[1])-m))
        print((b,A,n,diff),flush=True)
        assert diff<1e-11
        out.append({"base":b,"digits":A,"alpha":alpha,"c":c,"n":n,"depth":depth,
                    "explicit_gaps":gaps,"max_absolute_difference":diff})
    Path("independent-check.json").write_text(json.dumps({"mp_dps":mp.mp.dps,"rows":out},indent=2)+"\n")


if __name__=="__main__":
    main()
