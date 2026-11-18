"""High-precision fold roots and Gaussian saddle-weighted moment ratios."""
import json
from pathlib import Path

import mpmath as mp


def fold(gamma,b):
    r=gamma/b**2
    def h(x):
        v=1/(r*x*x)
        return v*mp.exp(v)*mp.e1(v)
    def derivative(x):
        a=r*x*x
        return (1-h(x))/a
    root=mp.findroot(lambda x:2*x*derivative(x)-1,(mp.mpf('.5'),mp.mpf('.6')))
    rec=-mp.log(4*(root-root*root*h(root)))
    return root,rec


def main():
    mp.mp.dps=54
    gamma=mp.mpf(2)
    raw=[]
    for b in (8,16,32,64,256,1024):
        x,log_ratio=fold(gamma,mp.mpf(b))
        r=gamma/b**2
        raw.append(dict(b=b,scaled_fold=str(x),log_ratio=str(log_ratio),
                        second_log_scaled=str((log_ratio+r/4)/r**2)))
    nd,w=mp.gauss_quadrature(80,'legendre')
    moments=[]
    for n in (128,512,2048,8192):
        bn=mp.sqrt(n)
        ln=max(-12,8-bn)
        up=mp.mpf(12)
        half=(up-ln)/2
        mid=(up+ln)/2
        norm=mp.exp(n/2)*mp.power(n,-mp.mpf(n)/2)*mp.power(2,mp.mpf(n-1)/2)*mp.gamma(mp.mpf(n+1)/2)
        def eval_at(z):
            _,lr=fold(gamma,bn+z)
            log_gaussian=n*mp.log1p(z/bn)-bn*z-z*z/2
            return mp.exp(log_gaussian+n*lr)
        val=half*mp.fsum(wj*eval_at(mid+half*zj) for zj,wj in zip(nd,w))/norm
        moments.append(dict(order=n,truncated_weighted_fold_ratio=str(val),
                            asymptotic_limit=str(mp.exp(-gamma/4))))
    Path('fold_precision.json').write_text(json.dumps(dict(dps=mp.mp.dps,mpmath=mp.__version__,
                       gamma=str(gamma),d="1",folds=raw,moments=moments),indent=2)+"\n")
    print(moments)


if __name__=='__main__':
    main()
