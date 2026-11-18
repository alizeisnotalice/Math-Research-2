"""Float64/MP periodic integrals and fixed-high-mode cancellation solve."""
import json
from pathlib import Path
import sys

import mpmath as mp
import numpy as np
import scipy
from scipy.optimize import root


def table(n):
    x=2*np.pi*np.arange(n)/n
    frequencies=np.arange(1,9)[:,None]
    return np.cos(frequencies*x),np.sin(frequencies*x)


def harmonics(c,s,basis,count=3):
    co,si=basis
    y=np.asarray(c)@co+np.asarray(s)@si
    weights=np.exp(-0.5*np.arange(1,count+1)[:,None]**2*y)
    z=weights.sum(axis=1)
    vals=(np.einsum('ij,ij->i',weights,co[:count])+1j*np.einsum('ij,ij->i',weights,si[:count]))/z
    return vals


def fit_low(c,s,basis):
    c,s=np.array(c),np.array(s)
    n=basis[0].shape[1]
    initial=np.r_[c[:3],s[:3]]
    def values(u):
        cv=np.r_[u[:3],c[3:]]
        sv=np.r_[u[3:],s[3:]]
        f=harmonics(cv,sv,basis)
        return np.r_[f.real,f.imag]
    attempt=root(values,initial,tol=1e-11)
    z=values(attempt.x)
    return {'success':bool(np.max(np.abs(z))<1e-10),'nodes':n,'low':attempt.x.tolist(),
            'residualInfinityNorm':float(np.max(np.abs(z)))}


def high_precision_reference(basis):
    c=[0,0,0,0.05,0.04,0,0,0]
    s=[0,0,0,0,0,0,0.02,0]
    guess=fit_low(c,s,basis)
    assert guess['success']
    with mp.workdps(86):
        def tables(n):
            angles=[2*mp.pi*i/n for i in range(n)]
            rows=[]
            for t in angles:
                cos=tuple(mp.cos(k*t) for k in range(1,5))
                sin=tuple(mp.sin(k*t) for k in range(1,5))
                h=mp.mpf('0.05')*cos[3]+mp.mpf('0.04')*mp.cos(5*t)+mp.mpf('0.02')*mp.sin(7*t)
                rows.append((cos,sin,h))
            return rows
        def evaluate(u,rows,count):
            x=[mp.mpf(0)]*count;y=[mp.mpf(0)]*count;den=[mp.mpf(0)]*count
            for cos,sin,h in rows:
                g=h+sum(u[k]*cos[k]+u[k+3]*sin[k] for k in range(3))
                for k in range(count):
                    q=mp.exp(-mp.mpf((k+1)**2)*g/2)
                    x[k]+=q*cos[k];y[k]+=q*sin[k];den[k]+=q
            return [v/d for v,d in zip(x,den)]+[v/d for v,d in zip(y,den)]
        a=tables(384); b=tables(768)
        sol=mp.findroot(lambda *u:evaluate(u,a,3),[mp.mpf(str(v)) for v in guess['low']],tol=mp.mpf('1e-78'))
        out=evaluate(sol,b,4)
        small=max(abs(out[i]) for i in (0,1,2,4,5,6))
        n4=mp.sqrt(out[3]**2+out[7]**2)
        assert small<mp.mpf('1e-75') and n4>mp.mpf('0.01')
        return {'highCosine':c[3:],'highSine':s[3:],'lowCoefficients':list(map(lambda x:mp.nstr(x,76),sol)),
                'maxCoordinateAtFirstThree':mp.nstr(small,30),'fourthHarmonicMagnitude':mp.nstr(n4,55)}


def main():
    src=Path(sys.argv[1]);out=Path(sys.argv[2])
    short=table(2048);long=table(4096)
    checks=[]
    if src.exists():
        for row in json.loads(src.read_text())['rows']:
            hit=row['best'];c,s=hit['cosine'],hit['sine']
            a=harmonics(c,s,short);b=harmonics(c,s,long)
            assert float(np.max(np.abs(a-b)))<1e-12
            rho=np.abs(b)
            entry={'family':row['family'],'amplitude':row['amplitude'],'metal':hit['metrics'],
                   'cpuMagnitudes':rho.tolist(),'maxAbsoluteDifference':float(np.max(np.abs(rho-hit['metrics'][:3])))}
            if row['family']=='corrected':entry['fixedHighModeFit']=fit_low(c,s,long)
            checks.append(entry)
    ref=high_precision_reference(long)
    print('fixed high mode first 3 residual',ref['maxCoordinateAtFirstThree'],'fourth magnitude',ref['fourthHarmonicMagnitude'],flush=True)
    if checks:print('maximum GPU discrepancy',max(x['maxAbsoluteDifference'] for x in checks),flush=True)
    out.write_text(json.dumps({'python':sys.version,'scipy':scipy.__version__,'numpy':np.__version__,'mpmath':mp.__version__,
           'selectedCandidates':checks,'fixedHighModeReference':ref},indent=2)+'\n')


if __name__=='__main__': main()
