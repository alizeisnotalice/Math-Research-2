"""Independent erfcx/Quadpack and high-precision checks of both Laplace scans."""
import json
import math
from pathlib import Path
import sys

import mpmath as mp
import scipy
from scipy.integrate import quad
from scipy.special import erfcx


SQ=math.sqrt(2*math.pi)


def params(hit,joint):
    p=hit['profile']
    b,t=p['b'],p['t']
    l=hit['lambda'];off=hit['offset'];L=math.log(l)
    logCap=math.log(2*b/t)-0.5*math.log(2*math.pi*t)-b*b/(2*t)
    Ccrit=b*b-2*logCap
    capSlope=math.exp(logCap)
    if joint:
        beta=4;share=p['share']
        C=Ccrit-2*math.log(share)
        c=2*math.exp(-b*b/2)/(SQ*(1-share)*capSlope)
        eta=8*b+off
        D=-2*(b*(b*b+3-C)+(1-share)/share*2/c)
    else:
        beta=p['beta'];share=1;C=Ccrit+off
        c=1;eta=D=0
    return b,t,l,off,L,beta,share,C,c,eta,D


def check64(hit,joint):
    b,t,l,off,L,beta,share,C,c,eta,D=params(hit,joint)
    K=beta*L+C
    def tail(y):
        z=y-b;x=z/l
        R=b*z*z/(l+z)-(beta*math.log1p(x)-K*x*(2+x))/(2*(1+x)**2)
        R+=(eta*(L+math.log1p(x))+D)/(2*l*(1+x)**3)
        return math.exp(-y*y/2)/SQ*math.expm1(R)
    dt=quad(tail,-8,8,epsabs=1e-13,epsrel=1e-12)[0]/math.erf(math.sqrt(32))
    a=(t*l-b)/math.sqrt(2*t);z=(t*l+b)/math.sqrt(2*t)
    jc=l*l*t*math.sqrt(2*math.pi*t)/(4*b)*(erfcx(a)-erfcx(z))
    def edge(v):
        u=2*v/(c*l*l)
        return math.exp(-v)*math.expm1((l-b)*u-u*u/2)
    de=quad(edge,0,64,epsabs=1e-13,epsrel=1e-12)[0] if joint else 0
    m=math.log1p(share*dt+(1-share)*de)-math.log(jc)
    if not joint:m+=(2-beta/2)*L-off/2
    return {'logMargin':m,'logTailCorrection':math.log1p(dt),'logCapCorrection':math.log(jc),
            'tailFirstTerm':b*(b*b+1+beta/2-K),
            'edgeSlope':c if joint else None}


def mp_check(hit,joint):
    with mp.workdps(68):
        p=hit['profile']
        b,t=map(lambda k:mp.mpf(str(p[k])),('b','t'))
        logCap=mp.log(2*b/t)-mp.log(2*mp.pi*t)/2-b*b/(2*t)
        Ccrit=b*b-2*logCap
        l,off=map(lambda k:mp.mpf(str(hit[k])),('lambda','offset'))
        L=mp.log(l)
        theta=mp.mpf(str(hit.get("theta",0)))
        capSlope=mp.exp(logCap)
        if joint:
            beta=4;share=mp.mpf(str(p['share']))
            C=Ccrit-2*mp.log(share)
            c=2*mp.exp(-b*b/2)/(mp.sqrt(2*mp.pi)*(1-share)*capSlope)
            eta=8*b+off
            D=-2*(b*(b*b+3-C)+(1-share)/share*2/c)
        else:
            beta=mp.mpf(str(p['beta']));share=1;C=Ccrit+off
            c=1;eta=D=0
        def tail(y):
            z=y-b;r=l+z
            gap=2*b/r-(beta*mp.log(r)+C)/r**2+(eta*mp.log(r)+D)/r**3+theta*mp.log(r)**2/r**4
            action=l*l*gap/2-b*l+b*z+(beta*L+C)/2
            return mp.exp(-y*y/2)*mp.expm1(action)/mp.sqrt(2*mp.pi)
        dt=mp.quad(tail,[-8,-4,0,4,8])/mp.erf(mp.sqrt(32))
        q0=(t*l-b)/mp.sqrt(2*t);q1=(t*l+b)/mp.sqrt(2*t)
        jc=l*l*t*mp.sqrt(2*mp.pi*t)/(4*b)*(mp.exp(q0*q0)*mp.erfc(q0)-mp.exp(q1*q1)*mp.erfc(q1))
        if joint:
            def edge(v):
                u=2*v/(c*l*l)
                return mp.exp(-v)*mp.expm1((l-b)*u-u*u/2)
            de=mp.quad(edge,[0,4,16,64])
        else:de=0
        ans=mp.log1p(share*dt+(1-share)*de)-mp.log(jc)
        if not joint:ans+=(2-beta/2)*L-off/2
        return ans


def main():
    work=Path(sys.argv[1]); output=Path(sys.argv[2]); report={}
    for name,is_joint in [('sync.json',False),('joint_admissible/joint_admissible.json',True)]:
        p=work/name
        if not p.exists():continue
        hits={}
        for label,c in json.loads(p.read_text())['hits'].items():
            d=check64(c,is_joint);m=mp_check(c,is_joint)
            assert abs(d['logMargin']-float(m))<3e-9,(name,label,d,m)
            hits[label]={'parameters':c,'reference':d,'mp68LogMargin':mp.nstr(m,48),
                         'absoluteMetalError':abs(d['logMargin']-c['result'][0])}
            print(name,label,'log margin',d['logMargin'],'GPU error',hits[label]['absoluteMetalError'],flush=True)
        report[name]=hits
    asympt=[]
    for lam in (256,4096,65536,1048576):
        for off in (-0.5,0,0.5):
            item={'profile':{'b':0.75,'t':0.5,'share':0.5},'lambda':lam,'offset':off}
            z=mp_check(item,True)
            asympt.append({'lambda':lam,'etaOffset':off,'scaledGap':float(z)*4*lam/math.log(lam)})
    second=[]
    for lam in (4096,1048576,100000000):
        for theta in (-20,-16,0):
            item={'profile':{'b':0.75,'t':0.5,'share':0.5},'lambda':lam,'offset':0,'theta':theta}
            z=mp_check(item,True)
            second.append({'lambda':lam,'theta':theta,'scaledGap':float(z)*2*lam**2/math.log(lam)**2,
                           'target':8+theta/2})
    print('largest lambda second order',second[-3:],flush=True)
    output.write_text(json.dumps({'secondOrderJointCheck':second,'jointAsymptoticCheck':asympt,'scipy':scipy.__version__,'mpmath':mp.__version__,
       'python':sys.version,'checked':report},indent=2)+'\n')


if __name__=='__main__':main()
