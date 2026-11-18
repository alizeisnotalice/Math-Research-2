"""Independent mpmath bisection of the core-touching witness width."""
import json
from pathlib import Path

import mpmath as mp
from validation import tail, boundary_integral


def evaluate(a,m,b,r,lam,width):
    left=max(r,lam-width/2)
    right=left+width
    z=mp.sqrt(1+a*lam*lam)
    j=mp.exp(lam**2/(2*z*z))/z*(tail(-r*z+lam/z)-tail(r*z+lam/z))
    power=(1-m)*lam*lam/2
    # Compute left contributions by independent endpoint quadrature.
    p_far=boundary_integral(lam+r,power)
    if left==r:
        outer_left=mp.mpf(0)
    else:
        q=lam-left
        outer_left=mp.exp(power-q*q/2)/mp.sqrt(2*mp.pi)*mp.quad(
            lambda t:mp.exp(-q*t-t*t/2), [0,left-r])
    outer_right=boundary_integral(right-lam,power)
    plateau=mp.exp((1-b)*lam*lam/2)*(1-tail(lam-left)-tail(right-lam))
    return j+p_far+outer_left+outer_right+plateau


def main():
    mp.mp.dps=65
    r=mp.mpf(1)
    a=mp.mpf(1)
    b=mp.mpf(2)
    rows=[]
    for n in [32,128,512,2048]:
        lam=mp.mpf(n)
        cases=[('fixed_beta', mp.mpf(1),None),('fixed_beta',mp.mpf(5)/2,None)]
        for tau in [mp.mpf(0),mp.mpf(2)]:
            cases.append(('touch_transition',2*r-2*(mp.log(lam)+tau)/lam,tau))
        for mode,beta,tau in cases:
            m=beta/lam
            touch=2*(lam-r)
            ftouch=evaluate(a,m,b,r,lam,touch)
            lo,hi=touch-8,touch+8
            for _ in range(52):
                mid=(lo+hi)/2
                if evaluate(a,m,b,r,lam,mid)>1:
                    lo=mid
                else:
                    hi=mid
            w=(lo+hi)/2
            assert abs(evaluate(a,m,b,r,lam,w)-1)<mp.mpf('1e-10')
            row=dict(mode=mode,lam=str(lam),beta=str(beta),width=str(w),
                     envelope_at_touch=str(ftouch),is_constrained=w>touch)
            if tau is not None:
                row.update(tau=str(tau),asymptotic_envelope_at_touch=str(mp.exp(tau-r*r/2)/mp.sqrt(2*mp.pi)))
            rows.append(row)
    Path('validation_corner.json').write_text(json.dumps(dict(precision_digits=mp.mp.dps,rows=rows),indent=2)+"\n")
    for r in rows[-4:]:
        print(r)


if __name__=='__main__':
    main()
