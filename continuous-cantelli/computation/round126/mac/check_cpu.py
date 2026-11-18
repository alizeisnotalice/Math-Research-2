"""Independent completed-square CDF, adaptive quadrature, and interval certificates."""
import json
import math
from pathlib import Path
import sys

import mpmath as mp
from scipy.integrate import quad
import scipy

mp.mp.dps=64
mp.iv.dps=50
STEPS=(.125,.25,.40,.60,.85,1.15,1.60,2.25)

def params(hit,ctx):
    A,m,d,wL=map(lambda x:ctx.mpf(str(x)),hit['geometry'])
    wR,sigma,delta,p=map(lambda x:ctx.mpf(str(x)),hit['shape'])
    return A,p,sigma,((m-d,wL,ctx.mpf(0)),(m+d,wR,delta))

def quadratic(A,sigma,well,c,a,ctx):
    center,w,offset=well
    root=ctx.sqrt(A)
    invT=ctx.exp(-A/2)
    x1=root*invT
    x0=(A+root*center)*invT
    B=x1*x1+(1+a/root)**2/(sigma*sigma)
    L=-x0*x1-a-a*a/root
    C=-(x0*x0+offset)/2-a*(center+(c+offset)/root)-(a*a/2)*(1+(2*root*center+c+offset)/A)
    return C,L,B,w

def moment_cdf(q,order):
    C,L,B,w=q
    mu=L/B
    k=mp.sqrt(B/2)
    I=mp.exp(C+L*L/(2*B))*mp.sqrt(mp.pi/(2*B))*(mp.erf(k*(w-mu))-mp.erf(k*(-w-mu)))
    ans=[I]
    qL=mp.exp(C-L*w-B*w*w/2)
    qR=mp.exp(C+L*w-B*w*w/2)
    for j in range(order):
        current=(L*ans[j]+(j*ans[j-1] if j else 0)-(w**j*qR-(-w)**j*qL))/B
        ans.append(current)
    return ans

def moment_adaptive(q,order):
    C,L,B,w=map(float,q)
    out=[]
    for j in range(order+1):
        val=quad(lambda y:y**j*math.exp(C+L*y-B*y*y/2),-w,w,epsabs=1e-12,epsrel=1e-12)[0]
        out.append(mp.mpf(val))
    return out

def moment_interval(q,order,panels=1024):
    C,L,B,w=q
    width=2*w
    step=width/panels
    acc=mp.iv.mpf(0)
    for i in range(panels):
        y=-w+(i+mp.iv.mpf('.5'))*step
        acc+=mp.iv.exp(C+L*y-B*y*y/2)
    m1=abs(L)+B*w
    upper=mp.iv.exp(C+abs(L)*w)*(m1*m1+B)
    err=(width**3*upper/(24*panels*panels)).b
    I=step*acc+mp.iv.mpf([-err,err])
    ans=[I]
    qL=mp.iv.exp(C-L*w-B*w*w/2)
    qR=mp.iv.exp(C+L*w-B*w*w/2)
    for j in range(order):
        current=(L*ans[j]+(j*ans[j-1] if j else 0)-(w**j*qR-(-w)**j*qL))/B
        ans.append(current)
    return ans

def evaluate(hit,moments,ctx,shifts):
    A,p,sigma,wells=params(hit,ctx)
    root=ctx.sqrt(A)
    T=[moments(quadratic(A,sigma,w,0,0,ctx),4) for w in wells]
    I=sum(m[0] for m in T)
    c=ctx.log(A)+2*ctx.log(I/(ctx.sqrt(2*ctx.pi)*p))
    totalD=0
    totalD2=0
    totalF=0
    depths=[]
    for (center,w,offset),(I0,I1,I2,I3,I4) in zip(wells,T):
        d0=-center-(c+offset)/root
        d2=-1/(sigma*sigma*root)
        totalD+=d0*I0-I1+d2*I2
        totalD2+=d0*d0*I0-2*d0*I1+(1+2*d0*d2)*I2-2*d2*I3+d2*d2*I4
        totalF+=(1+(2*root*center+c+offset)/A)*I0+2*I1/root+I2/(sigma*sigma*A)
        u=-root*sigma*sigma
        if ctx is mp.iv:
            assert u.b < -w.b
            u=-w
        else:
            u=max(-w,min(w,u))
        depths.append(1+(2*root*(center+u)+c+offset+(u/sigma)**2)/A)
    F=totalF/I
    d=totalD/I
    d2=totalD2/I
    baseline=p*F
    variance=p*d2+p*p*d*d/(1-p)
    rows=[]
    for step in shifts:
        a=ctx.mpf(str(step))
        qp,qm=[p*sum(moments(quadratic(A,sigma,w,c,sign*a,ctx),0)[0] for w in wells)/I for sign in (1,-1)]
        if qp>=1 or qm>=1:
            rows.append({'lambda':step,'qPlus':qp,'qMinus':qm})
        else:
            rem=1-p
            lh=baseline+rem/a**2*ctx.log(rem*rem/((1-qp)*(1-qm)))
            one_plus=baseline+2*rem/a**2*(-a*p*d/rem-ctx.log((1-qp)/rem))
            one_minus=baseline+2*rem/a**2*(a*p*d/rem-ctx.log((1-qm)/rem))
            rows.append({'lambda':step,'qPlus':qp,'qMinus':qm,'hellinger':lh,
                         'tiltedJensenPlus':one_plus,'tiltedJensenMinus':one_minus})
    return {'baseline':baseline,'variance':variance,'minima':depths,'shiftRows':rows,'c':c,'meanD':d}

def limiting_certificate():
    v=mp.iv
    p=v.mpf(57)/250
    left=v.mpf(51)/25  # 2.04
    right=v.mpf(103)/50 # 2.06
    h=v.mpf(5)/4
    cosh=lambda x:(v.exp(x)+v.exp(-x))/2
    variance=p*right**2
    cap0=p*cosh(right)
    cap1=p/2*(v.exp(right**2/2)+v.exp(-v.mpf('.5')-right))
    witness=p*(cosh(h*left)-1)-(v.exp(h*h/2)-1)
    assert variance<1 and cap0<1 and cap1<1 and witness>0
    return {'varianceUpper':str(variance),'massCapSmallShifts':str(cap0),'massCapLargeShifts':str(cap1),
            'hellingerStrictGapLower':str(witness)}

def finite_narrow_certificate(hit):
    ctx=mp.iv
    A,p,sigma,wells=params(hit,ctx)
    bounds=evaluate(hit,moment_interval,ctx,[1.25])
    row=bounds['shiftRows'][0]
    assert bounds['baseline']<1 and bounds['variance']<1 and row['hellinger']>1
    root=ctx.sqrt(A)
    c=bounds['c']
    I=[moment_interval(quadratic(A,sigma,w,0,0,ctx),0)[0] for w in wells]
    largest_weight=max((mass/sum(I)).b for mass in I)
    D=[]
    for center,w,offset in wells:
        y=ctx.mpf([-w.b,w.b])
        D.append(-center-y-(c+offset+(y/sigma)**2)/root)
    caps=[]
    for side,sign in ((0,1),(1,-1)):
        other=1-side
        a=(sign*D[side]).b
        b=(-sign*D[other]).a
        F_pos=bounds['minima'][side].a
        F_neg=bounds['minima'][other].a
        assert a>0 and b>0 and F_pos>0 and F_neg>0
        small=p*largest_weight*max(ctx.mpf(2), (ctx.exp(a)+ctx.exp(-b)).b)
        large=p*largest_weight*(ctx.exp(a*a/(2*F_pos))+ctx.exp(-b-F_neg/2))
        assert small<1 and large<1
        caps.append({'direction':sign,'smallShiftCap':str(small),'largeShiftCap':str(large)})
    return {'outwardInterval':recursive(bounds,str),'allRealIndividualShiftCaps':caps}

def recursive(obj,encode):
    if isinstance(obj,dict):
        return {k:recursive(v,encode) for k,v in obj.items()}
    if isinstance(obj,list):
        return [recursive(v,encode) for v in obj]
    if isinstance(obj,(float,int,str)):
        return obj
    return encode(obj)

def numeric_errors(left,right):
    if isinstance(left,dict):
        return max(numeric_errors(left[k],right[k]) for k in left if k in right)
    if isinstance(left,list):
        return max(numeric_errors(a,b) for a,b in zip(left,right))
    return abs(float(left)-float(right))

def main():
    root,out=map(Path,sys.argv[1:3])
    scan=json.loads(root.read_text())
    checked={}
    for k,hit in scan['hits'].items():
        a=evaluate(hit,moment_cdf,mp,STEPS)
        b=evaluate(hit,moment_adaptive,mp,STEPS)
        check=numeric_errors(a,b)
        assert check<2e-7,(k,check)
        verified={'parameters':hit,'cdf64':recursive(a,lambda x:mp.nstr(x,36)),'cdfVsAdaptiveMaxError':check}
        if k.startswith('best_strict'):
            interval=evaluate(hit,moment_interval,mp.iv,[hit['overlap'][0]])
            row=interval['shiftRows'][0]
            assert interval['baseline']<1 and interval['variance']<1 and min(t.a for t in interval['minima'])>0
            assert row['qPlus']<1 and row['qMinus']<1 and row['hellinger']>1
            verified['outwardInterval']=recursive(interval,str)
        checked[k]=verified
    witness={'geometry':[16384,0,2.05,.01],'shape':[.01,1,0,.228]}
    narrow=evaluate(witness,moment_cdf,mp,[1.25])
    out.write_text(json.dumps({'mpmath':mp.__version__,'scipy':scipy.__version__,'checks':checked,
        'limitingNarrowWindows':limiting_certificate(),'finiteNarrowWindowReplay':recursive(narrow,str),
        'finiteNarrowWindowCertificate':finite_narrow_certificate(witness)},indent=2)+'\n')
if __name__=='__main__':
    main()
