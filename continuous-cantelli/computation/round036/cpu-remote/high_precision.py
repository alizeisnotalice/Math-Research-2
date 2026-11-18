"""Independent high-precision packet contour diagnostic at the exact coefficients."""
import json
from pathlib import Path
import platform

import mpmath as mp


def main():
    mp.mp.dps=68
    h=mp.mpf('0.5')
    a=mp.mpf(3)/8
    g=mp.sqrt(3)/2
    b=3*mp.sqrt(3)/8
    w=-mp.mpf(1)/4+mp.j*g/2
    zstar=g+mp.j*h
    delta=mp.mpf(79)/100
    radius=12*mp.sqrt(3)
    q0=1+delta*zstar**4/2
    second=-3*zstar*zstar-delta/q0
    expected=mp.sqrt(2*mp.pi)/(mp.sqrt(q0)*mp.sqrt(-second))
    def phase(z):
        return mp.j*z-z**4/4+delta*(mp.j*z-z**4/2-w)**2/(2*(1+delta*z**4/2))
    def action(z):
        return z*z*(a+mp.re(phase(z)))
    peak=mp.findroot(lambda x:mp.diff(action,x),(mp.mpf('.7'),mp.mpf('.8')))
    rows=[]
    for lam in (32,64,128,256,512,1024):
        def integrand(x):
            z=x+mp.j*h
            return mp.exp(lam*(phase(z)+a-mp.j*b))/mp.sqrt(1+delta*z**4/2)
        core=mp.sqrt(lam)*mp.quad(integrand,[-3,-g,0,g-mp.mpf('.25'),g,g+mp.mpf('.25'),3])
        full=core+mp.sqrt(lam)*(mp.quad(integrand,[-radius,-4,-3])+
                               mp.quad(integrand,[3,4,radius]))
        rows.append(dict(lam=lam,scaled_packet=str(full),saddle_relative_error=str(abs(full/expected-1)),
                         outer_contour_difference=str(abs(full-core))))
    out=dict(python=platform.python_version(),mpmath=mp.__version__,digits=mp.mp.dps,
             delta=str(delta),real_action_stationary_point=str(peak),real_action_at_point=str(action(peak)),
             quadratures=rows)
    Path('high_precision.json').write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    main()
