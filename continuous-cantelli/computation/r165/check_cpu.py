"""Cantor geometry, independent normal-CDF replay, and certified center/split onset brackets."""
from collections import defaultdict
import json
import math
from pathlib import Path
import platform
import sys

import mpmath as mp
import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.special import ndtr
import scipy

mp.iv.dps=46
THETA_REF=math.sqrt(8/7)*math.exp(-2/7)
TARGETS=[j/100 for j in range(70,100)]+[THETA_REF]
SELECTED=(1/256,1/64,.125,.25,.5,1.,1.375)

def geometry(c,depth):
    length=1.
    centers=np.array([0.])
    for n in range(1,depth+1):
        next_length=length*(1-c*2.**(-n))/2
        shift=(length-next_length)/2
        centers=np.ravel(np.column_stack((centers-shift,centers+shift)))
        length=next_length
    tail=math.prod(1-c*2.**(-n) for n in range(depth+1,54))
    return centers,length,tail

def value(geom,v,x):
    centers,width,tail=geom
    sigma=math.sqrt(v)
    i=np.searchsorted(centers,x-9*sigma-width/2)
    j=np.searchsorted(centers,x+9*sigma+width/2)
    mid=(centers[i:j]-x)/sigma
    left,right=np.abs(mid)-width/(2*sigma),np.abs(mid)+width/(2*sigma)
    return tail*float(np.sum(ndtr(-left)-ndtr(-right)))

def curvature(geom,v):
    centers,width,tail=geom
    a,b=centers-width/2,centers+width/2
    return tail*float(np.sum(a*np.exp(-a*a/(2*v))-b*np.exp(-b*b/(2*v))))/(v*math.sqrt(2*math.pi*v))

def split_root(c,geom):
    va=np.geomspace(max((c/4)**2*.7,1e-8),.25,90)
    vv=[curvature(geom,v) for v in va]
    brackets=[(float(a),float(b)) for a,b,x,y in zip(va,va[1:],vv,vv[1:]) if x>0 and y<0]
    return brentq(lambda v:curvature(geom,v),*brackets[-1],xtol=5e-15)

def iv_geometry(c,N):
    ctx=mp.iv
    c=ctx.mpf(str(c))
    L=ctx.mpf(1)
    centers=[ctx.mpf(0)]
    for n in range(1,N+1):
        nxt=L*(1-c/(2**n))/2
        shift=(L-nxt)/2
        centers=[x for m in centers for x in (m-shift,m+shift)]
        L=nxt
    tail=ctx.mpf(1)
    M=N+60
    for n in range(N+1,M+1):
        tail*=1-c/(2**n)
    p=ctx.mpf([(tail*(1-c/(2**M))).a,tail.b])
    return centers[len(centers)//2:],L,p

def derivative_poly(n,order,v):
    poly={2*n:mp.iv.mpf(1)}
    for _ in range(order):
        nxt=defaultdict(lambda:mp.iv.mpf(0))
        for i,a in poly.items():
            if i:nxt[i-1]+=a*i
            nxt[i+1]-=a/v
        poly=dict(nxt)
    return poly

def eval_poly(p,x):
    return sum(a*x**j for j,a in p.items())

def derivative_bound(n,order,v):
    return sum(abs(a)*mp.iv.mpf('.5')**j for j,a in derivative_poly(n,order,v).items())

def moment_enclosure(geom,v,max_order):
    ctx=mp.iv
    v=v if hasattr(v,'_mpi_') else ctx.mpf(str(v))
    centers,L,tail=geom
    count=len(centers)
    mass=count*L
    weighted=[ctx.exp(-x*x/(2*v)) for x in centers]
    output=[]
    for n in range(max_order+1):
        poly2=derivative_poly(n,2,v)
        approximate=L*sum(g*(x**(2*n)+L*L/24*eval_poly(poly2,x)) for x,g in zip(centers,weighted))
        quaderr=(count*L**5*derivative_bound(n,4,v)/1920).b
        homogen=((1-tail)*mass*L*L*derivative_bound(n,2,v)/4).b
        err=(quaderr+homogen).b
        output.append(tail*approximate+ctx.mpf([-err,err]))
    return output

def certification(c,root):
    low=math.floor((root-1e-7)*1e9)/1e9
    high=math.ceil((root+1e-7)*1e9)/1e9
    g=iv_geometry(c,11)
    v0,v1=(mp.iv.mpf(str(v)) for v in (low,high))
    low_m=moment_enclosure(g,low,1)
    pos=low_m[1]-v0*low_m[0]
    assert pos>0,(c,low,pos)
    highest=math.ceil((.25/high-1)/2)
    ms=moment_enclosure(g,high,highest)
    coeff=[ms[n+1]-(2*n+1)*v1*ms[n] for n in range(highest)]
    assert all(x<0 for x in coeff),(c,high,coeff)
    assert mp.iv.mpf('.25')<=(2*highest+1)*v1
    height_hi=2*ms[0]/mp.iv.sqrt(2*mp.iv.pi*v1)
    whole=mp.iv.mpf([str(low),str(high)])
    whole_cut=math.ceil((.25/low-1)/2)
    both=moment_enclosure(g,whole,whole_cut)
    pitchfork=[both[n+1]-(2*n+1)*whole*both[n] for n in range(1,whole_cut)]
    assert all(x<0 for x in pitchfork),(c,whole,pitchfork)
    assert mp.iv.mpf('.25')<=(2*whole_cut+1)*v0
    height_at_root=2*both[0]/mp.iv.sqrt(2*mp.iv.pi*whole)
    return {'vSplitBelow':low,'vUnimodalAbove':high,'positiveMomentAtLower':str(pos),
            'negativeSuccessiveMomentDifferences':list(map(str,coeff)),
            'negativeHigherSuccessiveMomentsThroughoutBracket':list(map(str,pitchfork)),
            'centerHeightAtUpper':str(height_hi),'heightAtFirstSplit':str(height_at_root)}

def interval_reference_level(c,root):
    g=iv_geometry(c,10)
    theta=mp.iv.sqrt(mp.iv.mpf(8)/7)*mp.iv.exp(-mp.iv.mpf(2)/7)
    bracket=(math.floor((root-1e-7)*1e9)/1e9,math.ceil((root+1e-7)*1e9)/1e9)
    vals=[]
    for v in bracket:
        M0=moment_enclosure(g,v,0)[0]
        vals.append(2*M0/mp.iv.sqrt(2*mp.iv.pi*mp.iv.mpf(str(v)))-theta)
    assert vals[0]>0 and vals[1]<0
    return {'theta':str(theta),'vBracket':bracket,'centerHeightMinusTheta':list(map(str,vals))}

def global_center_cap_from_first_gap(c):
    a=c/4;b=.5
    v=(b*b-a*a)/(2*math.log(b/a))
    return 2*(ndtr(b/math.sqrt(v))-ndtr(a/math.sqrt(v)))

def strict_split_all_v_certificate(c):
    ctx=mp.iv
    a=ctx.mpf(str(c))/4;b=ctx.mpf('.5')
    v=(b*b-a*a)/(2*ctx.log(b/a))
    # Integrate over [a,b] at the exact peak variance with a midpoint enclosure.
    panels=2048
    dx=(b-a)/panels
    total=dx*sum(ctx.exp(-(a+dx*(j+ctx.mpf('.5')))**2/(2*v)) for j in range(panels))
    err=((b-a)**3*(1/v+b*b/(v*v))/(24*panels*panels)).b
    cap=2*(total+ctx.mpf([-err,err]))/ctx.sqrt(2*ctx.pi*v)
    assert cap<ctx.sqrt(ctx.mpf(8)/7)*ctx.exp(-ctx.mpf(2)/7)
    return {'c':c,'firstGapOptimizingVariance':str(v),'allVarianceCenterUpper':str(cap)}

def main():
    inpath,outpath,levelpath=map(Path,sys.argv[1:4])
    raw=json.loads(inpath.read_text()); grouped=defaultdict(list)
    for point in raw['rows']:grouped[point['c']].append(point)
    thresholds={}
    for c,rows in grouped.items():
        levels=[]
        for t in TARGETS:
            idx=max((i for i,row in enumerate(rows) if row['maxValue']>t),default=-1)
            if idx<0 or idx>=len(rows)-1:
                levels.append({'theta':t,'unbracketed':True})
            else:
                row=rows[idx]
                levels.append({'theta':t,'sampledBracketV':[row['v'],rows[idx+1]['v']],
                    'maxPointOnLowerVGrid':row['x'],'peakMinusCenterAtLowerVGrid':row['maxValue']-row['center'],
                    'varianceScaledCenterCurvatureAtLowerVGrid':row['centerCurvature']})
        thresholds[str(float(c))]=levels
    levelpath.write_text(json.dumps({'thetaReference':THETA_REF,'levels':thresholds},indent=2)+'\n')
    replay=[]
    for c in SELECTED:
        geom=geometry(c,14)
        fine=geometry(c,16)
        r=split_root(c,geom)
        on=certification(c,r)
        tests=[]
        for t in (THETA_REF,.99):
            g=next(q for q in thresholds[str(float(c))] if q['theta']==t)
            if 'sampledBracketV' not in g:continue
            v=g['sampledBracketV'][0]
            x=g['maxPointOnLowerVGrid']
            raw_row=next(q for q in grouped[c] if q['v']==v)
            vh,vc,xh=(value(fine,v,x),value(fine,v,0),x)
            match=abs(vh-value(geom,v,x))
            if x>0.004:
                rad=max(.004,math.sqrt(v)/3)
                best=minimize_scalar(lambda x:-value(fine,v,x),method='bounded',
                        bounds=(max(0,x-rad),min(.5,x+rad)),
                        options={'xatol':1e-12})
                xh=best.x
            tests.append({'theta':t,'gridV':v,'gpuPeak':raw_row['maxValue'],'gpuX':x,
                          'cpuAtGpuX':vh,'cpuAtZero':vc,'cpuLocalRefineX':xh,'cpuLocalRefineHeight':value(fine,v,xh),
                          'absoluteDifferenceStage14Vs16':match,'cpuMinusGpuAtGpuX':vh-raw_row['maxValue']})
        replay.append({'c':c,'splitRootApproximate':r,'heightAtSplitApproximate':value(fine,r,0),
                       'certifiedFirstSplit':on,'thresholdReplays':tests,
                       'firstGapOnlyCenterCap':global_center_cap_from_first_gap(c)})
    sample_ref=next(x for x in replay if x['c']==1/256)
    profile=geometry(1/256,16)
    ref=brentq(lambda v:value(profile,v,0)-THETA_REF,.05,.3,xtol=3e-15)
    reference=interval_reference_level(1/256,ref)
    assert reference['vBracket'][0]>sample_ref['certifiedFirstSplit']['vUnimodalAbove']
    outpath.write_text(json.dumps({'python':platform.python_version(),'scipy':scipy.__version__,'numpy':np.__version__,
        'mpmath':mp.__version__,'thetaReference':THETA_REF,'referenceCenterThreshold':reference,
        'cPoint25AllVarianceSplitAtReferenceThreshold':strict_split_all_v_certificate(.25),'checked':replay},indent=2)+'\n')

if __name__=='__main__':main()
