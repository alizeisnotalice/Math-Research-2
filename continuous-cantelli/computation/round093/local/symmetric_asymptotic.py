from pathlib import Path
import json

import mpmath as mp

mp.mp.dps=65


def profile(r,t):
    s=mp.sqrt(t)
    sf=lambda z:mp.erfc(z/mp.sqrt(2))/2
    phi=lambda z:mp.exp(-z*z/2)/mp.sqrt(2*mp.pi)
    p1=sf(r)
    base=2*r*(phi(r)-r*p1)
    centers=[(1,-4*r*n) for n in range(-3,4)]+[(-1,-(4*n+2)*r) for n in range(-3,4)]

    def moment(left,order):
        terms=[]
        for sign,c in centers:
            lo,hi=(left-c)/s,(r-c)/s
            mass=sf(lo)-sf(hi) if lo>=0 else sf(-hi)-sf(-lo)
            if order==0:
                terms.append(sign*mass)
                continue
            pl,ph=phi(lo),phi(hi)
            second=c*c*mass+2*c*s*(pl-ph)+t*(mass+lo*pl-hi*ph)
            terms.append(sign*second)
        return mp.fsum(terms)

    lo,hi=mp.mpf(0),mp.mpf(r)
    for _ in range(190):
        mid=(lo+hi)/2
        if moment(mid,0)>p1:
            lo=mid
        else:
            hi=mid
    cut=(lo+hi)/2
    gain=2*(r*r*p1-moment(cut,2))
    assert abs(moment(cut,0)/p1-1)<mp.mpf('1e-48')
    return dict(r=r,t=str(t),cut_over_r=mp.nstr(cut/r,22),
                allowed_mass=mp.nstr(2*p1,22),
                correction_over_base=mp.nstr(gain/base,22),
                asymptotic_ratio=mp.nstr(gain/(base*(1-t)*r*r),22))


def main():
    results=[profile(r,mp.mpf(t)) for t in ['.1','.25','.5','.9'] for r in [4,8,16,32,64]]
    Path(__file__).with_name('symmetric_asymptotic.json').write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps(results))


if __name__=='__main__':
    main()
