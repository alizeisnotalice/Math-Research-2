from math import ceil, log, sqrt
from pathlib import Path
import json

import numpy as np
from scipy.optimize import brentq
from scipy.special import iv

from multiscale_mellin import rule

ROOT=Path(__file__).parent


def point(l,beta,omega,a,d,panelscale):
    c=.5*log(2)-log(iv(0,a))
    initial=lambda b:l-log(2*b)

    def root(b):
        ell=initial(b)
        return b-(ell+.5*log(ell)+c+d*log(ell)/ell)

    b=brentq(root,l-2*log(l),l+2*log(l))
    ell_b=initial(b)
    freq=omega*beta*ell_b**(beta-1)/sqrt(2*b)
    prefactor=np.exp(c+d*log(ell_b)/ell_b)*sqrt(ell_b/(4*np.pi*b))
    w,weights=rule(-12,12,panelscale*ceil(24*(1.4+freq)),24)
    shift=-np.log1p(w/sqrt(2*b))
    lr=np.log1p(shift/ell_b)
    theta=(omega*ell_b**beta)%(2*np.pi)
    phase=theta+omega*ell_b**beta*np.expm1(beta*lr)
    ell=ell_b+shift
    delta=shift+.5*lr+d*(np.log(ell)/ell-log(ell_b)/ell_b)+a*np.cos(phase)
    penalty=(sqrt(2*b)*w-2*delta)**2/(4*(b+delta))
    return float(prefactor*(weights@np.exp(delta-penalty)))


def main():
    rows=[]
    for beta in [1.7,1.85]:
        for l in [4096,16384,65536,262144,1048576]:
            for d in [0,.25,.5]:
                f=point(l,beta,.9,.8,d,1)
                fine=point(l,beta,.9,.8,d,2)
                assert abs(f-fine)<1e-10,(l,beta,d,f,fine)
                rows.append(dict(L=l,beta=beta,d=d,value=f,
                                 scaled=l*(f-1)/log(l),
                                 residual=l*(f-1)-(d-.25)*log(l),
                                 grid_error=abs(f-fine)))
    (ROOT/'drift.json').write_text(json.dumps(rows,indent=2)+'\n')


if __name__=='__main__':
    main()
