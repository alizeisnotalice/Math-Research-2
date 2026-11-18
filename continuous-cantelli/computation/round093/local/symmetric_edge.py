from pathlib import Path
import json

import mpmath as mp

mp.mp.dps=62


def profile(r,t):
    var=mp.mpf(t)
    free_edge=mp.exp(-r*r/(2*var))/mp.sqrt(2*mp.pi*var)
    gaussian_tail=mp.erfc(r/mp.sqrt(2))/2
    target=r*gaussian_tail/free_edge
    scale=mp.sqrt(var*target)
    shifts=[4*n+1 for n in range(-3,4)]
    prefactors=[mp.exp(-(k*k-1)*r*r/(2*var)) for k in shifts]

    def p_scaled(v):
        u=scale*v
        paired=mp.fsum(c*mp.sinh(k*u/var) for c,k in zip(prefactors,shifts))
        return 2*mp.exp(-u*u/(2*var*r*r))*paired*scale/target

    v=mp.findroot(lambda end:mp.quad(p_scaled,[0,end])-1,(mp.mpf('.7'),mp.mpf('1.4')))
    assert abs(mp.quad(p_scaled,[0,v])-1)<mp.mpf('1e-54')
    m_corr=r*(mp.exp(-r*r/2)/mp.sqrt(2*mp.pi)/gaussian_tail-r)
    gain=scale*mp.quad(lambda z:(2*z-scale*z*z/r**2)*p_scaled(z),[0,v])/m_corr
    if var==1:
        predicted=6*mp.acosh(mp.mpf('1.5'))-2*mp.sqrt(5)
    else:
        predicted=mp.mpf(4)/3*var**mp.mpf('.75')*mp.exp(-(1-1/var)*r*r/4)
    return dict(r=r,t=str(var),rescaled_width=mp.nstr(scale*v,24),
                correction_over_base=mp.nstr(gain,24),predicted=mp.nstr(predicted,24),
                asymptotic_ratio=mp.nstr(gain/predicted,24))


def main():
    rows=[profile(r,t) for t in ['1','1.03','1.25','2','5'] for r in [4,8,16,32,64]]
    Path(__file__).with_name('symmetric_edge.json').write_text(json.dumps(rows,indent=2)+'\n')
    print(json.dumps(rows))


if __name__=='__main__':
    main()
