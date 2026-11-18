"""Independent adaptive/65-digit evaluation and outward certificate for a fixed bump pair."""
import json
import math
from pathlib import Path
import platform
import struct
import sys

import mpmath as mp
from scipy.integrate import quad
from scipy.optimize import brentq, minimize_scalar
import scipy

mp.mp.dps=65
mp.iv.dps=48

def f32(x):
    return struct.unpack('f',struct.pack('f',x))[0]

def profile_at(shape, origin, lam, side, high=False):
    ctx=mp if high else math
    ln=mp.mpf(lam) if high else float(lam)
    c,L,s=origin[:3]
    def integral(fn,lo,hi):
        if high:return mp.quad(fn,[lo,hi])
        return quad(fn,float(lo),float(hi),epsabs=1e-12,epsrel=1e-12)[0]
    def weight(u):
        if side=='mass':return 1
        a=u if side=='left' else 1-u
        return ctx.sin(ln*a) if ln else a
    def spatial(u):
        return c+L*(u-.5)
    denominator=integral(lambda u:weight(u)*ctx.exp(-spatial(u)**2/(2*s))/ctx.sqrt(s),0,1)
    value=0
    rat=mp.mpf('0.05') if high else .05
    for j in (0,1):
        wid,ht=shape[j],shape[j+2]
        def numerator(t):
            u=wid*t if j==0 else 1-wid*t
            rel=min(1,t/rat,(1-t)/rat)
            return wid*weight(u)*ctx.exp(-spatial(u)**2/2+ht*ln**2*rel)
        value+=sum(integral(numerator,a,b) for a,b in ((0,rat),(rat,1-rat),(1-rat,1)))
    return value/denominator

def exact_certificate():
    ctx=mp.iv
    k=ctx.mpf(14)/5
    osc_hi=ctx.exp(ctx.mpf(1)/8+ctx.mpf(3)/4*ctx.pi**2)*(1-ctx.cos(ctx.pi/100))
    osc_lo=ctx.exp(ctx.mpf(3)/4*k*k-ctx.mpf(1)/8)*(ctx.sin(k*999/2000)-ctx.sin(k*981/2000))/ctx.sin(k/2)
    mass=ctx.exp(ctx.mpf(1)/8)/50
    bump_mass=ctx.mpf(3)/(100*ctx.sqrt(2*ctx.pi))
    tent_at_point1=ctx.mpf(1)/10*ctx.exp(-ctx.mpf(1)/200)/ctx.sqrt(2*ctx.pi)
    assert osc_hi<1 and osc_lo>1 and mass<1 and bump_mass<tent_at_point1
    return {'endpointPiRatioUpper':str(osc_hi),'interiorK14over5RatioLower':str(osc_lo),
            'endpointZeroMassAndLinearUpper':str(mass),'addedGaussianMassUpper':str(bump_mass),
            'tentGaussianMassAtRadius1over10Lower':str(tent_at_point1)}

def main():
    inpath,outpath=map(Path,sys.argv[1:3])
    raw=json.loads(inpath.read_text()); rows=[]
    for bucket in raw['rows']:
        for tag in ('fixed','winner','selectedPhase'):
            if tag not in bucket:continue
            hit=bucket[tag]
            shape=list(map(f32,hit['shape']))
            origin=list(map(f32,hit['origin']))
            hshape,horigin=list(map(mp.mpf,shape)),list(map(mp.mpf,origin))
            cases=[('mass',0),('left',0),('right',0),('left',mp.pi),('right',mp.pi),
                ('left',f32(hit['sweep'][2])),('right',f32(hit['sweep'][3]))]
            results=[]
            for side,lam in cases:
                value=profile_at(shape,origin,lam,side)
                high=profile_at(hshape,horigin,lam,side,True)
                err=abs(value-float(high))
                assert err<1e-9,(bucket,tag,side,lam,err)
                results.append({'phase':side,'lambda':float(lam),'quadpack':value,'mp65':mp.nstr(high,45),'difference':err})
            opt=[]
            delta=math.pi/(raw['steps']-1)
            for j,side in enumerate(('left','right')):
                pivot=f32(hit['sweep'][j+2])
                polished=minimize_scalar(lambda v:-profile_at(shape,origin,v,side),method='bounded',
                    bounds=(max(0,pivot-delta),min(math.pi,pivot+delta)),options={'xatol':1e-12})
                opt.append({'phase':side,'lambda':polished.x,'ratio':-polished.fun})
            p=[x['quadpack'] for x in results]
            endpoints=[*p[:3],max(p[3:5])]
            diffs=[abs(f32(a)-b) for a,b in zip(hit['endpoint'],endpoints)]
            diffs.extend(abs(f32(a)-b) for a,b in zip(hit['sweep'][:2],p[5:7]))
            assert max(endpoints)<1 and max(p[5:7])>1,(bucket,tag,p)
            rows.append({'mode':bucket['mode'],'kind':tag,'input':hit,'results':results,'polished':opt,
                         'maxAbsoluteGpuCpuDifference':max(diffs)})
    shape=[mp.mpf(1)/100]*2+[mp.mpf(3)/4]*2
    origin=[mp.mpf(0),mp.mpf(1),mp.mpf(1)]
    ratios={name:mp.nstr(profile_at(shape,origin,lam,side,True),46)
        for name,lam,side in (('mass',0,'mass'),('linear',0,'left'),('pi',mp.pi,'left'),('interior',mp.mpf(14)/5,'left'))}
    def peak_height(x):
        if abs(x)<.49 or abs(x)>.5:return 0
        t=(.5-abs(x))/.01
        return 1.5*min(1,t/.05,(1-t)/.05)
    B=2*sum(quad(lambda x:peak_height(x)*math.exp(-x*x/2)/math.sqrt(2*math.pi),a,b)[0]
            for a,b in ((.49,.4905),(.4905,.4995),(.4995,.5)))
    def deficit(radius):
        return 2*quad(lambda x:(1-x/radius)*math.exp(-x*x/2)/math.sqrt(2*math.pi),0,radius)[0]
    radius=brentq(lambda r:deficit(r)-B,.001,.1,xtol=1e-15)
    outpath.write_text(json.dumps({'python':platform.python_version(),'scipy':scipy.__version__,'mpmath':mp.__version__,
        'exactWitness':exact_certificate(),'witnessRatios':ratios,'tentRadiusApproximate':radius,'checked':rows},indent=2)+'\n')
if __name__=='__main__':main()
