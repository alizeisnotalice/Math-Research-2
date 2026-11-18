"""High-precision numerical values of the proven all-cut first-order boundary."""
import json
import platform
import mpmath as mp

mp.mp.dps=55
rows=[]
for H in map(mp.mpf, ('1.001','1.01','1.1','1.5','2','4','8')):
    def difference(q):
        tau=1-mp.gammainc(mp.mpf('1.5'),0,q)/(q*mp.gammainc(mp.mpf('.5'),0,q))
        return H*mp.expm1(q)/mp.expm1(H*q)-tau
    lo,hi=mp.mpf('1e-25'),mp.mpf('1')
    while difference(hi)>0:
        hi*=2
    for _ in range(170):
        mid=(lo+hi)/2
        if difference(mid)>0:
            lo=mid
        else:
            hi=mid
    root=(lo+hi)/2
    easy=mp.log(mp.mpf('1.5'))/(H+1)
    assert root>easy
    rows.append(dict(H=mp.nstr(H,6),all_cut_range_upper=mp.nstr(root,22),simple_bound=mp.nstr(easy,14)))
print(json.dumps(dict(python=platform.python_version(),mpmath=mp.__version__,decimal_precision=mp.mp.dps,
                      bisection_iterations=170,seed=None,rows=rows),indent=2))
