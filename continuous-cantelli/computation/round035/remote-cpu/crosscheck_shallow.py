"""Independent mpmath dual integral and two-dimensional Gaussian quadrature."""
import json
from pathlib import Path

import mpmath as mp
from crosscheck import gaussian_dual, phi


def main():
    mp.mp.dps=48
    a,w,p,h,offset,half=map(mp.mpf,['8','0.4420026','1','2.9063306','.25','.03125'])
    lower,upper=a+offset-half,a+offset+half
    exact_dual=gaussian_dual(lower,upper,[a-w,a,a+w],
                            lambda x:1-h/a**2*(1-((x-a)/w)**2)**p)
    nodes,weights=mp.gauss_quadrature(40,'legendre')
    primal=mp.fsum(weights[i]*weights[j]*phi(w*nodes[i]-offset-half*nodes[j])
                   *mp.exp(h*(1-nodes[i]**2)**p*(1+(offset+half*nodes[j])/a)**2/2)
                   for i in range(len(nodes)) for j in range(len(nodes)))*w/2
    assert abs(primal/exact_dual-1)<mp.mpf('1e-40')
    out=dict(a=str(a),w=str(w),p=str(p),h=str(h),offset=str(offset),halfwidth=str(half),
             precision=mp.mp.dps,primal=str(primal),dual=str(exact_dual),
             relative_discrepancy=str(abs(primal/exact_dual-1)))
    Path('crosscheck_shallow.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    main()
