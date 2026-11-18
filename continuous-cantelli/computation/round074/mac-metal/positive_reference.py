"""High-precision direct positive kernel on selected Cantor-gap profile."""
import json
from pathlib import Path
import sys
import mpmath as mp


def main():
    original=json.loads(Path(sys.argv[1]).read_text())
    row=min(original['rows'],key=lambda item:item['best6']['scores'][0])
    hit=row['best6']
    B=row['base'];D=[j for j in range(B) if (hit['mask']>>j)&1];depth=hit['depth']
    gaps=[(u+1,v) for u,v in zip(D,D[1:]) if v>u+1]
    values=[]
    with mp.workdps(58):
        alpha,h=mp.mpf(str(row['exponent'])),mp.mpf(str(hit['height']))
        density=mp.mpf(len(D))/B
        for n in range(1,7):
            arg=2*mp.pi*n
            gain=mp.mpf(n*n)*h/2
            phase=mp.mpc(1);live=mp.mpf(1);scale=mp.mpf(1)
            C=mp.mpc(0);Z=mp.mpf(0)
            for _ in range(depth):
                for lo,hi in gaps:
                    g=mp.mpf(hi-lo)*scale/B
                    center=mp.mpf(lo+hi)*scale/(2*B)
                    def integrand(t,oscillation):
                        z=t**4/2
                        off=g*(z-mp.mpf('0.5'))
                        return 4*t**3*mp.exp(-gain*(g*z)**alpha)*(mp.cos(arg*off) if oscillation else 1)
                    V=mp.quad(lambda t:integrand(t,True),[0,0.25,0.5,1])*live*g/scale
                    W=mp.quad(lambda t:integrand(t,False),[0,0.25,0.5,1])*live*g/scale
                    C+=phase*mp.exp(1j*arg*center)*V
                    Z+=W
                phase*=sum(mp.exp(1j*arg*mp.mpf(d)*scale/B) for d in D)/len(D)
                live*=density;scale/=B
            C+=phase*live*mp.exp(1j*arg*scale/2)*mp.sinc(arg*scale/2)
            Z+=live
            tail=live*gain*(scale/2)**alpha
            values.append({'n':n,'real':mp.nstr(C.real,45),'imag':mp.nstr(C.imag,35),
                           'magnitude':mp.nstr(abs(C),45),'positiveWeightMass':mp.nstr(Z,45),
                           'normalizedByPositiveWeight':mp.nstr(abs(C)/Z,45),
                           'normalizedByLossWeight':mp.nstr(abs(C)/(1-Z),45),
                           'finiteToInfiniteAbsoluteBound':mp.nstr(tail,35)})
    out={'base':B,'digits':D,'depth':depth,'exponent':row['exponent'],'height':hit['height'],
         'mpmath':mp.__version__,'precisionDigits':58,'harmonics':values}
    Path(sys.argv[2]).write_text(json.dumps(out,indent=2)+'\n')
    print('base',B,'digits',D,'depth',depth)
    for row in values: print(row['n'],row['magnitude'],row['normalizedByPositiveWeight'],row['finiteToInfiniteAbsoluteBound'],flush=True)


if __name__=='__main__':main()
